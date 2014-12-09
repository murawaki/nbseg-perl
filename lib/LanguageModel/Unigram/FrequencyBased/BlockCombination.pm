package LanguageModel::Unigram::FrequencyBased::BlockCombination::CombinationStruct;
#
# data structure for combination of frequency counts
#
use strict;
use warnings;
use utf8;

use Math::GSL::RNG qw/gsl_rng_uniform_int/;

use LanguageModel::Util qw/$rng randPartition/;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	total => 0,
	freq2idList => {}, # freq -> [id] : holders sorted by frequency
	# freqList : freq -> list of assignment combinations
    };
    bless($self, $class);
    return $self;    
}

sub append {
    my ($self, $freq, $id) = @_;
    $self->{total} += $freq;
    {
	no warnings qw/uninitialized/;
	push(@{$self->{freq2idList}->{$freq}}, $id);
    }
}

sub freqList {
    my ($self) = @_;
    return $self->{freqList};
}

sub setFreqList {
    my ($self) = @_;

    my $candList = { 0 => [[]] };
    while ((my ($freq, $idList) = each(%{$self->{freq2idList}}))) {
	my $candList2 = {};
	while ((my ($totalFreq, $statusList) = each(%$candList))) {
	    foreach my $status (@$statusList) {
		push(@{$candList2->{$totalFreq}}, $status);
		foreach my $i (1 .. scalar(@$idList)) {
		    my $clone = [];
		    map { push(@$clone, $_); } (@$status);
		    push(@$clone, [$freq, $i]);
		    push(@{$candList2->{$totalFreq + $freq * $i}}, $clone);
		}
	    }
	}
	$candList = $candList2;
    }
    $self->{freqList} = $candList;
}

sub randSelect {
    my ($self, $freq) = @_;

    my $list = $self->{freqList}->{$freq};
    return [] unless(scalar(@$list) > 0);

    my $asgnCmb;
    if (scalar(@$list) == 1) {
	$asgnCmb = $list->[0];
    } else {
	# select the assignment combination
	my $massList = [];
	foreach my $l (@$list) {
	    my $mass = 1;
	    foreach my $tmp (@$l) {
		my ($eFreq, $count) = @$tmp;
		my $total = scalar(@{$self->{freq2idList}->{$eFreq}});
		$mass *= &combination($total, $count);
	    }
	    push(@$massList, $mass);
	}
	my $selected = randPartition($massList);
	$asgnCmb = $list->[$selected];
    }

    my $rv = [];
    foreach my $l (@$asgnCmb) {
	my ($eFreq, $count) = @$l;
	my $idList = $self->{freq2idList}->{$eFreq};
	my $total = scalar(@$idList);

	my $clone = [];
	map { push(@$clone, $_) } @$idList;
	while (scalar(@$clone) > $total - $count) {
	    my $r = gsl_rng_uniform_int($rng->raw, scalar(@$clone));
	    my $selected = splice(@$clone, $r, 1);
	    push(@$rv, $selected);
	}
    }
    return $rv;
}

sub total {
    my ($self) = @_;
    return $self->{total};
}

sub combination {
    my ($n, $r) = @_;
    my $product = 1;
    while ($r > 0) {
        $product *= $n--;
        $product /= $r--;
    }
    return $product;
}

1;

package LanguageModel::Unigram::FrequencyBased::BlockCombination;
#
# type-based sampler for text in which each sentence has a frequency count
# type-based sampling over the site across multiple sentences
#
use strict;
use warnings;
use utf8;

use Carp::Assert;
use Math::Cephes qw/lgam/;
use Scalar::Util qw/refaddr/;

use LanguageModel::Util qw/$rng randPartitionLog/;
use Text::Sentence;  # cast compounds as sentences

use base qw/LanguageModel::Unigram::TypeBased/;

# our $MAXLGM = 2.556348e305; # DEBUG
our ($ALPHA_A, $ALPHA_B) = (0.01, 0.01); # G(a,b), the prior for the hyperparameter

#
# frozen text only increases word counts
# does not modify $self->{lexicon}->{$name}->{list}
#

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	zerogram => shift,
	opt => shift,
	unigram => {},      # word := count
	word2posList => {}, # word -> sid -> wpos := [sentence, node]
	tokenCount => 0,
	skipBlocks => {},   # signature 'left+right'
	totalRuns => 0,
	acceptedCount => 0, # Gibbs equivalent of Metropolis-Hastings
	alpha => 100000,
	a => $ALPHA_A, b => $ALPHA_B,
    };
    $self->{alpha} = $self->{opt}->{alpha} if (defined($self->{opt}->{alpha}));
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
    my $freq = $sentence->freq;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $word = $sentence->word($node);
	$self->addWord($word, $freq);
	$self->addPosList($word, $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $freq = $sentence->freq;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $word = $sentence->word($node);
	$self->removeWord($word, $freq);
	$self->removePosList($word, $sid, $iter->wpos) unless ($isFrozen);
    }
}

sub sampleSentence {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    my $sid = $sentence->sid;
    my $freq = $sentence->freq;
    my $iter = $sentence->charIterator;
    while ((my $region = $iter->next)) {
	my $blockName = $region->{left} . '+' . $region->{right};
	next if (defined($self->{skipBlocks}->{$blockName}));
	$self->{skipBlocks}->{$blockName} = 1;
	next if ($self->doSkip($region));

	my $segType = $self->sampleBoundary($region, $iter, $anneal, $freq);
    }
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal, $freq) = @_;

    my ($leftCount, $rightCount, $segmentedList) = $self->getSegmentedList($region, $iter);
    my ($unsegmentedCount, $unsegmentedList) = $self->getUnsegmentedList($region, $iter);
    my $numList = scalar(@$segmentedList) + scalar(@$unsegmentedList);

    # use simpler sampling procedure if no other site is found
    return $self->sampleOneBoundary($region, $iter, $anneal) if ($numList <= 1);
    $self->{totalRuns}++;

    my $llength = length($region->{left});
    foreach my $holder (@$segmentedList) {
	my $sid = $holder->{sentence}->sid;
	my $lwpos = $holder->wpos;
	my $freq = $holder->{sentence}->freq;
	$self->removeWord($region->{left}, $freq);
	$self->removeWord($region->{right}, $freq);
	$self->removePosList($region->{left}, $sid, $lwpos);
	$self->removePosList($region->{right}, $sid, $lwpos + $llength);
    }
    foreach my $holder (@$unsegmentedList) {
	my $sid = $holder->{sentence}->sid;
	my $wpos = $holder->wpos;
	my $freq = $holder->{sentence}->freq;
	$self->removeWord($region->{unsegmented}, $freq);
	$self->removePosList($region->{unsegmented}, $sid, $wpos);
    }

    my $cmbStruct = $self->expandCandList($segmentedList, $unsegmentedList);
    my $candList = $cmbStruct->freqList;
    my $blockTotal = $cmbStruct->total;
    my @freqList = sort { $a <=> $b } (keys(%$candList));

    my $ALPHA = $self->{alpha};
    my $unsegmentedCharProb = $self->zeroProb($region->{unsegmented});
    my $leftCharProb = $self->zeroProb($region->{left});
    my $rightCharProb = $self->zeroProb($region->{right});

    my $isSymmetric = ($region->{left} eq $region->{right})? 1 : 0;
    my $logProbList = {};

    # $i: the number of unsegmented regions
    foreach my $i (@freqList) {
	my $logProb;
	if (defined($logProbList->{$i})) {
	    # use symmetry to reduce computation
	    $logProb = $logProbList->{$i};
	} else {
	    # sCm == s! / (m! * (s-m)!); ignore invariant s!
	    $logProb = -1 * (lgam($i + 1) + lgam($blockTotal - $i + 1));
	    $logProbList->{$i} = $logProb;
	    if (defined($candList->{$blockTotal - $i})
		&& !defined($logProbList->{$blockTotal - $i})) {
		$logProbList->{$blockTotal - $i} = $logProb;
	    }
	}

	# totalCount
	$logProb = -1 * lgam($ALPHA + $self->{tokenCount} + $blockTotal * 2 - $i);
	$logProb += lgam($ALPHA * $unsegmentedCharProb + ($unsegmentedCount + $i));
	if ($isSymmetric) {
	    # NOTE the calculation order; the first term might be by far smallest
	    $logProb += lgam($ALPHA * $leftCharProb + ($leftCount + 2 * ($blockTotal - $i)));
	} else {
	    $logProb += lgam($ALPHA * $leftCharProb + ($leftCount + $blockTotal - $i));
	    $logProb += lgam($ALPHA * $rightCharProb + ($rightCount + $blockTotal - $i));
	}
	$logProbList->{$i} += $logProb;	
    }

    my $status; # ordered list of region to be unsegmented
    {
	my $logProbTable;
	foreach my $freq (@freqList) {
	    push(@$logProbTable, $logProbList->{$freq});
	}
	my $idx = randPartitionLog($logProbTable, $anneal);
	my $selected = $cmbStruct->randSelect($freqList[$idx]);
	my @tmp = sort { $a <=> $b } (@$selected);
	$status = \@tmp;
    }

    # now update the real sentences
    my $newSegType;

    my $current = shift(@$status); $current = $blockTotal unless (defined($current));
    for (my $i = 0; $i < $numList; $i++) {
	my $type = Text::Sentence->BOUNDARY;
	if ($i == $current) {
	    $type = Text::Sentence->NONE_BOUNDARY;
	    $current = shift(@$status); $current = $blockTotal unless (defined($current));
	}

	if ($i < scalar(@$segmentedList)) {
	    my $holder = $segmentedList->[$i];
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    my $freq = $holder->{sentence}->freq;
	    if ($holder->{node} == $region->{node}) {
		$newSegType = $type;
	    }
	    if ($type == Text::Sentence->BOUNDARY) {
		$self->addWord($region->{left}, $freq);
		$self->addWord($region->{right}, $freq);
		$self->addPosList($region->{left}, $sid, $wpos, $holder->{sentence}, $holder->{node});
		$self->addPosList($region->{right}, $sid, $wpos + $llength, $holder->{sentence}, $holder->{node}->[Text::Sentence->NEXT]);
	    } else {
		my $newNode = $holder->merge($region->{unsegmented});
		$self->addWord($region->{unsegmented}, $freq);
		$self->addPosList($region->{unsegmented}, $sid, $wpos, $holder->{sentence}, $newNode);
	    }
	} else {
	    my $holder = $unsegmentedList->[$i - scalar(@$segmentedList)];
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    my $freq = $holder->{sentence}->freq;
	    if ($holder->{node} == $region->{node}) {
		$newSegType = $type;
	    }
	    if ($type == Text::Sentence->BOUNDARY) {
		my $newNode = $holder->split($region->{left}, $region->{right});
		$self->addWord($region->{left}, $freq);
		$self->addWord($region->{right}, $freq);
		$self->addPosList($region->{left}, $sid, $wpos, $holder->{sentence}, $newNode);
		$self->addPosList($region->{right}, $sid, $wpos + $llength, $holder->{sentence}, $newNode->[Text::Sentence->NEXT]);
	    } else {
		$self->addWord($region->{unsegmented}, $freq);
		$self->addPosList($region->{unsegmented}, $sid, $wpos, $holder->{sentence}, $holder->{node});
	    }
	}
    }
    return $newSegType;
}

sub sampleOneBoundary {
    my ($self, $region, $iter) = @_;

    my $sentence = $iter->{sentence};
    my $sid = $sentence->sid;
    my $wpos = $iter->wpos;
    my $freq = $iter->{sentence}->freq;
    my $llength = length($region->{left});
    if ($region->{segType} == Text::Sentence->BOUNDARY) {
	$self->removePosList($region->{left}, $sid, $wpos, $sentence, $iter->{node});
	$self->removePosList($region->{right}, $sid, $wpos + $llength, $sentence, $iter->{node}->[Text::Sentence->NEXT]);
	$self->removeWord($region->{left}, $freq);
	$self->removeWord($region->{right}, $freq);
    } else {
	$self->removePosList($region->{unsegmented}, $sid, $wpos, $sentence, $iter->{node});
	$self->removeWord($region->{unsegmented}, $freq);
    }
    # frequency counts
    my $leftCount = $self->getCountAndPosList($region->{left});
    my $rightCount = $self->getCountAndPosList($region->{right});
    my $unsegmentedCount = $self->getCountAndPosList($region->{unsegmented});

    my $ALPHA = $self->{alpha};
    my $isSymmetric = ($region->{left} eq $region->{right})? 1 : 0;
    my $logProbUnseg = -1 * lgam($ALPHA + $self->{tokenCount} + $freq);
    my $logProbSeg = -1 * lgam($ALPHA + $self->{tokenCount} + $freq * 2);
    $logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{unsegmented}) + $unsegmentedCount + $freq);
    $logProbSeg += lgam($ALPHA * $self->zeroProb($region->{unsegmented}) + $unsegmentedCount);
    if ($isSymmetric) {
	$logProbSeg += lgam($ALPHA * $self->zeroProb($region->{left}) + $leftCount + 2 * $freq);
	$logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{left}) + $leftCount);
    } else {
	$logProbSeg += lgam($ALPHA * $self->zeroProb($region->{left}) + $leftCount + $freq);
	$logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{left}) + $leftCount);
	$logProbSeg += lgam($ALPHA * $self->zeroProb($region->{right}) + $rightCount + $freq);
	$logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{right}) + $rightCount);
    }

    if (randPartitionLog([$logProbSeg, $logProbUnseg]) == 0) {
	my $newNode = $iter->{node};
	if ($region->{segType} != Text::Sentence->BOUNDARY) {
	    $newNode = $iter->split($region->{left}, $region->{right});
	}
	$self->addPosList($region->{left}, $sid, $wpos, $sentence, $newNode);
	$self->addPosList($region->{right}, $sid, $wpos + $llength, $sentence, $newNode->[Text::Sentence->NEXT]);
	$self->addWord($region->{left}, $freq);
	$self->addWord($region->{right}, $freq);
	return Text::Sentence->BOUNDARY;
    } else {
	my $newNode = $iter->{node};
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->addPosList($region->{unsegmented}, $sid, $wpos, $sentence, $newNode);
	$self->addWord($region->{unsegmented}, $freq);
	return Text::Sentence->NONE_BOUNDARY;
    }
}

sub clearIteration {
    my ($self) = @_;

    $self->{skipBlocks} = {};

    if ($self->{opt}->{debug}) {
	printf STDERR ("accept ratio: %f (%d / %d)\n",
		       $self->{acceptedCount} / $self->{totalRuns},
		       $self->{acceptedCount}, $self->{totalRuns});
    }

    $self->{totalRuns} = 0;
    $self->{acceptedCount} = 0;
}

sub expandCandList {
    my ($self, $segmentedList, $unsegmentedList) = @_;

    # printf STDERR ("%s\n", 2**(scalar(@$segmentedList)+scalar(@$unsegmentedList)));
    my $cmbStruct = LanguageModel::Unigram::FrequencyBased::BlockCombination::CombinationStruct->new;
    my $i = 0;
    foreach my $list (($segmentedList, $unsegmentedList)) {
	foreach my $holder (@$list) {
	    my $freq = $holder->{sentence}->freq;
	    $cmbStruct->append($freq, $i);
	    $i++;
	}
    }
    $cmbStruct->setFreqList;
    return $cmbStruct;
}

# output: freq. count after removal, list
sub getUnsegmentedList {
    my ($self, $region, $iter) = @_;
    my ($unsegCount, $unsegList) = $self->getCountAndPosList($region->{unsegmented});
    my $list = [];
    my $count = 0;
    while ((my ($sid, $uList) = each(%$unsegList))) {
	while ((my ($wpos, $tmp) = each(%$uList))) {
	    my ($sentence, $node) = @$tmp;
	    if ($node == $region->{node}) {
		unshift(@$list, $iter);
	    } else {
		push(@$list, $sentence->getNodeHolder($node, $wpos));
	    }
	    $count -= $sentence->freq;
	}
    }
    return ($unsegCount, $list);
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
		    my $freq = $sentence->freq;
		    $leftCount -= $freq;
		    $rightCount -= $freq;

		    # ASSERT
		    unless (defined($lnode->[Text::Sentence->NEXT])) {
			use Dumpvalue;
			printf STDERR ("%s\n", $region->{left});
			printf STDERR ("%s\n", $sentence->asText);
			Dumpvalue->new->dumpValue($lnode);
		    }
		    unless ($lnode->[Text::Sentence->NEXT] == $rlist->{$lwpos + $l}->[1]) {
			printf STDERR ("right node is not placed next to left node: %s-%s (%s)\n",
				       $region->{left}, $region->{right}, $sentence->asText)
			    if ($self->{opt}->{debug});
		    }
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
	    $leftCount -= $iter->{sentence}->freq;
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
		$leftCount -= $sentence->freq;
		$registered->{refaddr($lnode)} = 1;
		$registered->{refaddr($rnode)} = 1;
	    }
	}
	$rightCount = $leftCount;  # for completeness
    }
    return ($leftCount, $rightCount, $list);
}

sub hasMutableLatent { return 0 }
sub sampleLatent {}
# succeeding UnigramModel#sampleHyper

1;
