package LanguageModel::Unigram::TypeBased;

use strict;
use warnings;
use utf8;

use Carp::Assert;
use Scalar::Util qw/refaddr/;
use Math::Cephes qw/lgam/;
use Math::GSL::Randist qw/gsl_ran_flat/;

use LanguageModel::Util qw/$rng randPartitionLog/;
use Text::Sentence;

use base qw /LanguageModel::Unigram/;

#
# frozen text only increases word counts
# does not modify $self->{lexicon}->{$name}->{list}
#

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	alpha => 0.1,
	a => 1, b => 1, # G(a,b), the prior for the hyperparameter
	zerogram => shift,
	opt => shift,
	unigram => {},      # word := count
	word2posList => {}, # word -> sid -> wpos := [sentence, node]
	tokenCount => 0,
	skipBlocks => {},   # signature 'left+right'
	totalRuns => 0,
	acceptedCount => 0, # Gibbs equivalent of Metropolis-Hastings
    };
    $self->{alpha} = $self->{opt}->{alpha} if (defined($self->{opt}->{alpha}));
    $self->{lambda} = $self->{opt}->{lambda} if (defined($self->{opt}->{lambda}));
    $self->{opt}->{matchingSkip} = 0 unless (defined($self->{opt}->{matchingSkip}));
    $self->{opt}->{doSampleHyper} = 0 unless (defined($self->{opt}->{doSampleHyper}));
    bless($self, $class);
    return $self;
}

# sid must be assigned in advance
sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $word = $sentence->word($node);
	$self->addWord($word);
	$self->addPosList($word, $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $word = $sentence->word($node);
	$self->removeWord($word);
	$self->removePosList($word, $sid, $iter->wpos) unless ($isFrozen);
    }
}

sub getCountAndPosList {
    my ($self, $name) = @_;
    my $count = $self->{unigram}->{$name} || 0;
    my $posList = $self->{word2posList}->{$name} || {};
    return wantarray ? ($count, $posList) : $count;
}

sub addPosList {
    my ($self, $name, $sid, $wpos, $sentence, $node) = @_;

    no warnings qw/uninitialized/;
    assert(!defined($self->{word2posList}->{$name}->{$sid}->{$wpos}),
	   sprintf("override existing node: %s %d-%d\n", $name, $sid, $wpos)) if (DEBUG);
    $self->{word2posList}->{$name}->{$sid}->{$wpos} = [$sentence, $node];
}

sub removePosList {
    my ($self, $name, $sid, $wpos) = @_;

    my $posList = $self->{word2posList}->{$name};
    assert(defined($posList),
	   sprintf("%s not found in position list (%d-%d)\n", $name, $sid, $wpos)) if (DEBUG);
    assert(defined($posList->{$sid}->{$wpos}),
	   sprintf("%s not found in position %d-%d\n", $name, $sid, $wpos)) if (DEBUG);
    delete($posList->{$sid}->{$wpos});
    delete($posList->{$sid}) if (scalar(keys(%{$posList->{$sid}})) <= 0);
    # NOTE: we may have frozen count; count is greater than or equal to the number of positions
    if (scalar(keys(%$posList)) <= 0) {
	delete($self->{word2posList}->{$name});
    }
}

sub sampleSentence {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    my $iter = $sentence->charIterator;
    while ((my $region = $iter->next)) {
	my $blockName = $region->{left} . '+' . $region->{right};
	next if (defined($self->{skipBlocks}->{$blockName}));
	$self->{skipBlocks}->{$blockName} = 1;
	next if ($self->doSkip($region));

	my $segType = $self->sampleBoundary($region, $iter, $anneal);
    }
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    my ($leftCount, $rightCount, $segmentedList) = $self->getSegmentedList($region, $iter);
    my ($unsegmentedCount, $unsegmentedList) = $self->getUnsegmentedList($region, $iter);
    my $numState = scalar(@$segmentedList) + scalar(@$unsegmentedList);

    # use simpler sampling procedure if no other site is found
    return $self->sampleOneBoundary($region, $iter, $anneal) if ($numState <= 1);
    $self->{totalRuns}++;

    my $llength = length($region->{left});
    foreach my $holder (@$segmentedList) {
	my $sid = $holder->{sentence}->sid;
	my $lwpos = $holder->wpos;
	$self->removeWord($region->{left});
	$self->removeWord($region->{right});
	$self->removePosList($region->{left}, $sid, $lwpos);
	$self->removePosList($region->{right}, $sid, $lwpos + $llength);
    }
    foreach my $holder (@$unsegmentedList) {
	my $sid = $holder->{sentence}->sid;
	my $wpos = $holder->wpos;
	$self->removeWord($region->{unsegmented});
	$self->removePosList($region->{unsegmented}, $sid, $wpos);
    }

    my $ALPHA = $self->{alpha};
    my $logProbTable = [];
    # $i: the number of unsegmented regions
    for (my $i = 0, my $l = int($numState / 2 + 0.5); $i <= $l; $i++) {
	# sCm == s! / (m! * (s-m)!); ignore invariant s!
	my $logProb = -1 * (lgam($i + 1) + lgam($numState - $i + 1));
	$logProbTable->[$i] = $logProbTable->[$numState - $i] = $logProb;
    }

    my $unsegmentedCharProb = $self->zeroProb($region->{unsegmented});
    my $leftCharProb = $self->zeroProb($region->{left});
    my $rightCharProb = $self->zeroProb($region->{right});

    my $isSymmetric = ($region->{left} eq $region->{right})? 1 : 0;
    for (my $i = 0; $i <= $numState; $i++) {
	# totalCount
	my $logProb = -1 * lgam($ALPHA + $self->{tokenCount} + $numState * 2 - $i);
	$logProb += lgam($ALPHA * $unsegmentedCharProb + ($unsegmentedCount + $i));
	if ($isSymmetric) {
	    # NOTE the calculation order; the first term might be by far smallest
	    $logProb += lgam($ALPHA * $leftCharProb + ($leftCount + 2 * ($numState - $i)));
	} else {
	    $logProb += lgam($ALPHA * $leftCharProb + ($leftCount + $numState - $i));
	    $logProb += lgam($ALPHA * $rightCharProb + ($rightCount + $numState - $i));
	}
	$logProbTable->[$i] += $logProb;
    }
    my $m = randPartitionLog($logProbTable, $anneal);
    if (scalar(@$unsegmentedList) != $m) {
	$self->{acceptedCount}++; # changed
    }

    my $rest = $m;
    my $newSegType;
    for (my $i = 0; $i < $numState; $i++) {
	my $type = Text::Sentence->BOUNDARY; # segment
	if (gsl_ran_flat($rng->raw, 0, $numState - $i) < $rest) {
	    $type = Text::Sentence->NONE_BOUNDARY;
	    $rest--;
	}
	if ($i < scalar(@$segmentedList)) {
	    my $holder = $segmentedList->[$i];
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    if ($holder->{node} == $region->{node}) {
		$newSegType = $type;
	    }
	    if ($type == Text::Sentence->BOUNDARY) {
		$self->addWord($region->{left});
		$self->addWord($region->{right});
		$self->addPosList($region->{left}, $sid, $wpos, $holder->{sentence}, $holder->{node});
		$self->addPosList($region->{right}, $sid, $wpos + $llength, $holder->{sentence}, $holder->{node}->[Text::Sentence->NEXT]);
	    } else {
		my $newNode = $holder->merge($region->{unsegmented});
		$self->addWord($region->{unsegmented});
		$self->addPosList($region->{unsegmented}, $sid, $wpos, $holder->{sentence}, $newNode);
	    }
	} else {
	    my $holder = $unsegmentedList->[$i - scalar(@$segmentedList)];
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    if ($holder->{node} == $region->{node}) {
		$newSegType = $type;
	    }
	    if ($type == Text::Sentence->BOUNDARY) {
		my $newNode = $holder->split($region->{left}, $region->{right});
		$self->addWord($region->{left});
		$self->addWord($region->{right});
		$self->addPosList($region->{left}, $sid, $wpos, $holder->{sentence}, $newNode);
		$self->addPosList($region->{right}, $sid, $wpos + $llength, $holder->{sentence}, $newNode->[Text::Sentence->NEXT]);
	    } else {
		$self->addWord($region->{unsegmented});
		$self->addPosList($region->{unsegmented}, $sid, $wpos, $holder->{sentence}, $holder->{node});
	    }
	}
    }
    return $newSegType;
}

sub sampleOneBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    my $sentence = $iter->{sentence};
    my $sid = $sentence->sid;
    my $wpos = $iter->wpos;
    my $llength = length($region->{left});
    if ($region->{segType} == Text::Sentence->BOUNDARY) {
	$self->removePosList($region->{left}, $sid, $wpos, $sentence, $iter->{node});
	$self->removePosList($region->{right}, $sid, $wpos + $llength, $sentence, $iter->{node}->[Text::Sentence->NEXT]);
    } else {
	$self->removePosList($region->{unsegmented}, $sid, $wpos, $sentence, $iter->{node});
    }
    my $newSegType = $self->SUPER::sampleBoundary($region, $iter, $anneal);
    if ($newSegType == Text::Sentence->BOUNDARY) {
	$self->addPosList($region->{left}, $sid, $wpos, $sentence, $iter->{node});
	$self->addPosList($region->{right}, $sid, $wpos + $llength, $sentence, $iter->{node}->[Text::Sentence->NEXT]);
    } else {
	$self->addPosList($region->{unsegmented}, $sid, $wpos, $sentence, $iter->{node});
    }
    return $newSegType;
}

sub clearIteration {
    my ($self) = @_;
    $self->SUPER::clearIteration;

    $self->{skipBlocks} = {};

    # printf STDERR ("accept ratio: %f (%d / %d)\n",
    # 		   $self->{acceptedCount} / $self->{totalRuns},
    # 		   $self->{acceptedCount}, $self->{totalRuns});

    $self->{totalRuns} = 0;
    $self->{acceptedCount} = 0;
}

# output: freq. count after removal, list
sub getUnsegmentedList {
    my ($self, $region, $iter) = @_;
    my ($unsegCount, $unsegList) = $self->getCountAndPosList($region->{unsegmented});
    my $list = [];
    while ((my ($sid, $uList) = each(%$unsegList))) {
	while ((my ($wpos, $tmp) = each(%$uList))) {
	    my ($sentence, $node) = @$tmp;
	    if ($node == $region->{node}) {
		unshift(@$list, $iter);
	    } else {
		push(@$list, $sentence->getNodeHolder($node, $wpos));
	    }
	}
    }
    return ($unsegCount - scalar(@$list), $list);
}

sub getSegmentedList {
    my ($self, $region, $iter) = @_;

    # Following [Liang+ 2010], we first collect the left nodes and the right nodes separately,
    #   and then get the intersection
    my ($leftCount, $leftList) = $self->getCountAndPosList($region->{left});
    my ($rightCount, $rightList) = $self->getCountAndPosList($region->{right});
    my $list = [];
    my $l = length($region->{left});
    if ($region->{left} ne $region->{right}) {
	while ((my ($sid, $llist) = each(%$leftList))) {
	    my $rlist = $rightList->{$sid};
	    next unless (defined($rlist));
	    while ((my ($lwpos, $lstruct) = each(%$llist))) {
		my ($sentence, $lnode) = @$lstruct;
		next if ($sentence->segType($lnode) == Text::Sentence->FIXED_BOUNDARY);
		if (defined($rlist->{$lwpos + $l})) {
		    if ($lnode == $region->{node}) {
			unshift(@$list, $iter);
		    } else {
			push(@$list, $sentence->getNodeHolder($lnode, $lwpos));
		    }

		    assert(defined($lnode->[Text::Sentence->NEXT]),
			   sprintf("no right node for %s <- %s\n",
				   $region->{left},
				   $sentence->asText)) if (DEBUG);
		    assert($lnode->[Text::Sentence->NEXT] == $rlist->{$lwpos + $l}->[1],
			   sprintf("right node is not placed next to left node: %s-%s (%s)\n",
				   $region->{left}, $region->{right},
				   $sentence->asText)) if (DEBUG);
		}
	    }
	}
    } else {
	# need to check possible conflicts
	my $registered = {};
	my $shared = $leftList;  # NOTE: $leftList == $rightList
	if ($region->{segType} == Text::Sentence->BOUNDARY) {
	    # target nodes must always be sampled, so register them first
	    push(@$list, $iter);
	    $registered->{refaddr($iter->{node})} = 1;
	    $registered->{refaddr($iter->{node}->[Text::Sentence->NEXT])} = 1;
	}
	while ((my ($sid, $sharedList) = each(%$shared))) {
	    while ((my ($lwpos, $lstruct) = each(%$sharedList))) {
		my ($sentence, $lnode) = @$lstruct;
		next if ($sentence->segType($lnode) == Text::Sentence->FIXED_BOUNDARY);
		next if (defined($registered->{refaddr($lnode)}));
		my $rstruct = $sharedList->{$lwpos + $l};
		next unless (defined($rstruct));
		my (undef,     $rnode) = @$rstruct;
		next if (defined($registered->{refaddr($rnode)}));
		push(@$list, $sentence->getNodeHolder($lnode, $lwpos));
		$registered->{refaddr($lnode)} = 1;
		$registered->{refaddr($rnode)} = 1;
	    }
	}
    }
    return ($leftCount - scalar(@$list), $rightCount - scalar(@$list), $list);
}

sub hasMutableLatent { return 0 }
sub sampleLatent {}
# succeeding UnigramModel#sampleHyper

1;
