package LanguageModel::Bigram::TypeBased;

use strict;
use warnings;
use utf8;

use Carp::Assert;
use Scalar::Util qw/refaddr/;

use Text::Sentence;
use LanguageModel::Util qw/$rng shuffle randPartitionLog getFlipProposal/;
use LanguageModel::Table::DirichletProcess; # explicit table tracking for hyperparameter estimation

use base qw/LanguageModel::Bigram LanguageModel::Unigram::TypeBased/;

our $PROB_PERMUTATION = 0.1;
our $PROPOSAL_LAMBDA = 0.5;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	alpha0 => 10000,
	alpha1 => 10000,
	zerogramModel => shift,
	word2posList => {}, # word -> sid -> wpos := [sentence, node]
	opt => shift,
    };
    $self->{opt}->{type} = 'Dirichlet' unless (defined($self->{opt}->{type}));
    $self->{opt}->{matchingSkip} = 0 unless (defined($self->{opt}->{matchingSkip}));
    $self->{opt}->{doSampleHyper} = 0 unless (defined($self->{opt}->{doSampleHyper}));
    $self->{opt}->{nested} = 0 unless (defined($self->{opt}->{nested}));
    bless($self, $class);
    $self->init;
    return $self;
}

sub init {
    my ($self) = @_;

    my $tableClass = ($self->{opt}->{type} eq 'Dirichlet')?
	'LanguageModel::Table::DirichletProcess' : 'LanguageModel::Table::PitmanYorProcess';
    $self->{zerogram} = LanguageModel::Table::ZeroAdaptor->new($self->{zerogramModel}, !($self->{opt}->{nested}) );
    $self->{unigram} = $tableClass->new($self->{zerogram}, \&LanguageModel::Bigram::reduceUnigram);
    $self->{bigram} = $tableClass->new($self->{unigram}, \&LanguageModel::Bigram::reduceBigram);
    if ($self->{opt}->{type} eq 'Dirichlet') {
	$self->{unigram}->alpha($self->{opt}->{alpha} || $self->{alpha0});
	$self->{bigram}->alpha($self->{opt}->{alpha1} || $self->{alpha1});
    } else {
	$self->{unigram}->theta($self->{opt}->{alpha} || $self->{alpha0});
	$self->{bigram}->theta($self->{opt}->{alpha1} || $self->{alpha1});
    }
}

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->add($current, $prev);
	$self->addPosList($current, $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
	$prev = $current;
    }
    $self->{bigram}->add('$', $prev);
}

sub removeSentence {
    my ($self, $sentence, $docID, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->remove($current, $prev);
	$self->removePosList($current, $sid, $iter->wpos) unless ($isFrozen);
	$prev = $current;
    }
    $self->{bigram}->remove('$', $prev);
}

sub addPosList {
    my ($self, $name, $sid, $wpos, $sentence, $node) = @_;
    $self->LanguageModel::Unigram::TypeBased::addPosList($name, $sid, $wpos, $sentence, $node);
}

sub removePosList {
    my ($self, $name, $sid, $wpos, $sentence, $node) = @_;
    $self->LanguageModel::Unigram::TypeBased::removePosList($name, $sid, $wpos, $sentence, $node);
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

    $self->{totalRuns}++;

    # Implementation Note:
    #  We block-sample a unigram-level type block,
    #  but we should be aware of possible bigram-level conflicts
    #  Consider the sequence 'C A B A B D'.
    #  When sampling the block 'A B', the bigram-level regions are
    #  'C A B A' and 'B A B D' and they conflict with each other.
    #  As a workaround, we only sample one of these regions.
    my $arrangement = {};
    {
	# always sample the current region
	my $outLeftNode = $iter->{node}->[Text::Sentence->PREV];
	$arrangement->{refaddr($outLeftNode)} = 1 if (defined($outLeftNode));
	my $outRightNode = ($region->{segType} == Text::Sentence->NONE_BOUNDARY)?
	    $iter->{node}->[Text::Sentence->NEXT]
	    : $iter->{node}->[Text::Sentence->NEXT]->[Text::Sentence->NEXT];
	$arrangement->{refaddr($outRightNode)} = 1 if (defined($outRightNode));	
    }
    my ($leftCount, $rightCount, $segmentedList) = $self->getSegmentedList($region, $iter, $arrangement);
    my ($unsegmentedCount, $unsegmentedList) = $self->getUnsegmentedList($region, $iter, $arrangement);
    my $numState = scalar(@$segmentedList) + scalar(@$unsegmentedList);

    return $self->sampleOneBoundary($region, $iter, $anneal) if ($numState <= 1);

    my ($proposal, $logQF, $logQB) = getFlipProposal(scalar(@$unsegmentedList), $numState);
    my ($splitList, $mergeList) = $self->setArrangements($proposal, $numState, $segmentedList, $unsegmentedList);

    my $llength = length($region->{left});
    my $currentLogProb = 0.0;
    foreach my $holder (@$splitList) {
	# remove unseg counts
	my $current = $holder->{node};
	my $prevNode = $current->[Text::Sentence->PREV];
	my $nextNode = $current->[Text::Sentence->NEXT];
	my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	$self->{bigram}->remove($outRight, $region->{unsegmented});
	my $prob = $self->{bigram}->prob($outRight, $region->{unsegmented});
	$self->{bigram}->remove($region->{unsegmented}, $outLeft);
	$prob *= $self->{bigram}->prob($region->{unsegmented}, $outLeft);
	$currentLogProb += log($prob);
    }
    foreach my $holder (@$mergeList) {
	# remove left and right counts
	my $left = $holder->{node};
	my $prevNode = $left->[Text::Sentence->PREV];
	my $right = $left->[Text::Sentence->NEXT];
	my $nextNode = $right->[Text::Sentence->NEXT];
	my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	$self->{bigram}->remove($outRight, $region->{right});
	my $prob = $self->{bigram}->prob($outRight, $region->{right});
	$self->{bigram}->remove($region->{right}, $region->{left});
	$prob *= $self->{bigram}->prob($region->{right}, $region->{left});
	$self->{bigram}->remove($region->{left}, $outLeft);
	$prob *= $self->{bigram}->prob($region->{left}, $outLeft);
	$currentLogProb += log($prob);
    }

    my $proposalLogProb = 0.0;
    #  add bigram and unigram counts according to $splitList and $mergeList
    foreach my $holder (reverse(@$splitList)) {
	# add left and right counts
	my $current = $holder->{node}; # unsegmented
	my $prevNode = $current->[Text::Sentence->PREV];
	my $nextNode = $current->[Text::Sentence->NEXT];
	my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	my $prob = $self->{bigram}->prob($region->{left}, $outLeft);
	$self->{bigram}->add($region->{left}, $outLeft);
	$prob *= $self->{bigram}->prob($region->{right}, $region->{left});
	$self->{bigram}->add($region->{right}, $region->{left});
	$prob *= $self->{bigram}->prob($outRight, $region->{right});
	$self->{bigram}->add($outRight, $region->{right});
	$proposalLogProb += log($prob);
    }
    foreach my $holder (reverse(@$mergeList)) {
	# remove left and right counts
	my $left = $holder->{node};
	my $prevNode = $left->[Text::Sentence->PREV];
	my $right = $left->[Text::Sentence->NEXT];
	my $nextNode = $right->[Text::Sentence->NEXT];
	my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	my $prob = $self->{bigram}->prob($region->{unsegmented}, $outLeft);
	$self->{bigram}->add($region->{unsegmented}, $outLeft);
	$prob *= $self->{bigram}->prob($outRight, $region->{unsegmented});
	$self->{bigram}->add($outRight, $region->{unsegmented});
	$proposalLogProb += log($prob);
    }
    if (randPartitionLog([$proposalLogProb + $logQB, $currentLogProb + $logQF], $anneal) == 0) {
	# accepted
	# no need to change counts
	# modify split/merge nodes in sentences and add them to POS list
	$self->{acceptedCount}++;

	my $newSegType; # now check the state of the current node
	foreach my $holder (@$splitList) {
	    # add left and right counts
	    my $current = $holder->{node};
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    $self->removePosList($region->{unsegmented}, $sid, $wpos);
	    if ($current == $region->{node}) {
		$newSegType = Text::Sentence->BOUNDARY;
	    }
	    my $newLeft = $holder->split($region->{left}, $region->{right});
	    $self->addPosList($region->{left}, $sid, $wpos, $holder->{sentence}, $newLeft);
	    $self->addPosList($region->{right}, $sid, $wpos + $llength, $holder->{sentence}, $newLeft->[Text::Sentence->NEXT]);
	}
	foreach my $holder (@$mergeList) {
	    # remove left and right counts
	    my $left = $holder->{node};
	    my $sid = $holder->{sentence}->sid;
	    my $wpos = $holder->wpos;
	    $self->removePosList($region->{left}, $sid, $wpos);
	    $self->removePosList($region->{right}, $sid, $wpos + $llength);
	    if ($left == $region->{node}) {
		$newSegType = Text::Sentence->NONE_BOUNDARY;
	    }
	    my $newUnsegmented = $holder->merge($region->{unsegmented});
	    $self->addPosList($region->{unsegmented}, $sid, $wpos, $holder->{sentence}, $newUnsegmented);
	}
	return (defined($newSegType))? $newSegType : $region->{segType}; # might have been unaffected
    } else {
	# rejected
	# reverting add/remove operations
	# POS list unchanged
	foreach my $holder (@$splitList) {
	    # add left and right counts
	    my $current = $holder->{node};
	    my $prevNode = $current->[Text::Sentence->PREV];
	    my $nextNode = $current->[Text::Sentence->NEXT];
	    my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	    my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	    $self->{bigram}->remove($outRight, $region->{right});
	    $self->{bigram}->remove($region->{right}, $region->{left});
	    $self->{bigram}->remove($region->{left}, $outLeft);
	    $self->{bigram}->add($region->{unsegmented}, $outLeft);
	    $self->{bigram}->add($outRight, $region->{unsegmented});
	}
	foreach my $holder (@$mergeList) {
	    # remove left and right counts
	    my $left = $holder->{node};
	    my $prevNode = $left->[Text::Sentence->PREV];
	    my $right = $left->[Text::Sentence->NEXT];
	    my $nextNode = $right->[Text::Sentence->NEXT];
	    my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	    my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	    $self->{bigram}->remove($outRight, $region->{unsegmented});
	    $self->{bigram}->remove($region->{unsegmented}, $outLeft);
	    $self->{bigram}->add($region->{left}, $outLeft);
	    $self->{bigram}->add($region->{right}, $region->{left});
	    $self->{bigram}->add($outRight, $region->{right});
	}
	return $region->{segType}; # no change
    }
}

sub sampleOneBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    my $sentence = $iter->{sentence};
    my $sid = $sentence->sid;
    my $wpos = $iter->wpos;
    my $llength = length($region->{left});
    if ($region->{segType} == Text::Sentence->BOUNDARY) {
	$self->removePosList($region->{left}, $sid, $wpos);
	$self->removePosList($region->{right}, $sid, $wpos + $llength);
    } else {
	$self->removePosList($region->{unsegmented}, $sid, $wpos);
    }
    my $newSegType = $self->LanguageModel::Bigram::sampleBoundary($region, $iter, $anneal);
    if ($newSegType == Text::Sentence->BOUNDARY) {
	$self->addPosList($region->{left}, $sid, $wpos, $sentence, $iter->{node});
	$self->addPosList($region->{right}, $sid, $wpos + $llength, $sentence, $iter->{node}->[Text::Sentence->NEXT]);
    } else {
	$self->addPosList($region->{unsegmented}, $sid, $wpos, $sentence, $iter->{node});
    }
    return $newSegType;
}

sub getUnsegmentedList {
    my ($self, $region, $iter, $arrangement) = @_;
    my ($unsegCount, $listTmp) = $self->LanguageModel::Unigram::TypeBased::getUnsegmentedList($region, $iter);

    # now checks bigram-level conflicts
    my $list = [];
    if ($region->{segType} == Text::Sentence->NONE_BOUNDARY) {
	assert($iter == $listTmp->[0], sprintf("something wrong with unseg list: %s\n", $region->{unsegmented})) if (DEBUG);
	shift(@$listTmp);
	push(@$list, $iter);
    }
    foreach my $holder (@$listTmp) {
	my $current = $holder->{node};
	if (defined($arrangement->{refaddr($current)})) {
	    $unsegCount++;
	} else {
	    push(@$list, $holder);
	    my $outLeftNode = $holder->{node}->[Text::Sentence->PREV];
	    my $outRightNode = $holder->{node}->[Text::Sentence->NEXT];
	    $arrangement->{refaddr($outLeftNode)} = 1 if (defined($outLeftNode));
	    $arrangement->{refaddr($outRightNode)} = 1 if (defined($outRightNode));
	}
    }
    return ($unsegCount, $list);
}

sub getSegmentedList {
    my ($self, $region, $iter, $arrangement) = @_;
    my ($leftCount, $rightCount, $listTmp) = $self->LanguageModel::Unigram::TypeBased::getSegmentedList($region, $iter);

    # now checks bigram-level conflicts
    my $list = [];
    if ($region->{segType} == Text::Sentence->BOUNDARY) {
	assert($iter == $listTmp->[0], sprintf("something wrong with unseg list: %s\n", $region->{unsegmented})) if (DEBUG);
	shift(@$listTmp);
	push(@$list, $iter);
    }
    foreach my $holder (@$listTmp) {
	my $left = $holder->{node};
	my $right = $left->[Text::Sentence->NEXT];
	if (defined($arrangement->{refaddr($left)})
	    || defined($arrangement->{refaddr($right)})) {
	    $leftCount++;
	    $rightCount++;
	} else {
	    push(@$list, $holder);
	    my $outLeftNode = $holder->{node}->[Text::Sentence->PREV];
	    my $outRightNode = $holder->{node}->[Text::Sentence->NEXT]->[Text::Sentence->NEXT];
	    $arrangement->{refaddr($outLeftNode)} = 1 if (defined($outLeftNode));
	    $arrangement->{refaddr($outRightNode)} = 1 if (defined($outRightNode));
	}
    }
    return ($leftCount, $rightCount, $list);
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

sub setArrangements {
    my ($self, $proposal, $numState, $segmentedList, $unsegmentedList) = @_;

    my $array = [];
    foreach (1 .. $numState - $proposal) { push(@$array, 0); }
    foreach (1 .. $proposal) { push(@$array, 1); }
    if ($proposal != scalar(@$unsegmentedList)) {
	shuffle($array);
    } else {
	do {
	    shuffle($array);
	} while ($self->isIdentical($array, $numState - $proposal));
    }

    my ($splitList, $mergeList) = ([], []);
    foreach my $i (0 .. scalar(@$segmentedList) - 1) {
	if(shift(@$array) == 1) {
	    push(@$mergeList, $segmentedList->[$i])
	}
    }
    foreach my $i (0 .. scalar(@$unsegmentedList) - 1) {
	if(shift(@$array) == 0) {
	    push(@$splitList, $unsegmentedList->[$i])
	}
    }
    return ($splitList, $mergeList);
}

sub isIdentical {
    my ($self, $array, $zeroSize) = @_;
    for my $i (0 .. scalar(@$array) - 1) {
	if ($i < $zeroSize) {
	    return 0 if($array->[$i] != 0);
	} else {
	    return 0 if($array->[$i] != 1);
	}
    }
    return 1;
}

1;
