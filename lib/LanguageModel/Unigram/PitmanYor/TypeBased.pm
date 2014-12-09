package LanguageModel::Unigram::PitmanYor::TypeBased;
#
# the unigram language model based on a Pitman-Yor process
#   with hybrid type-based sampling
#
use strict;
use warnings;
use utf8;

use Text::Sentence;
use LanguageModel::Unigram;
use LanguageModel::Util qw/shuffle randPartitionLog getFlipProposal/;
use LanguageModel::Table::PitmanYorProcess;

use base qw/LanguageModel::Unigram::PitmanYor LanguageModel::Unigram::TypeBased/;

our $NULL_CHAR = '__NULL__';

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	zerogramModel => shift,
	word2posList => {}, # word -> sid -> wpos := [sentence, node]
	skipBlocks => {},   # signature 'left+right'
	opt => shift,
    };
    $self->{opt}->{doSampleHyper} = 0 unless (defined($self->{opt}->{doSampleHyper}));
    $self->{opt}->{nested} = 0 unless (defined($self->{opt}->{nested}));
    bless($self, $class);
    $self->init;
    return $self;
}

sub init {
    my ($self) = @_;

    $self->{zerogram} = LanguageModel::Table::ZeroAdaptor->new($self->{zerogramModel}, !($self->{opt}->{nested}) );
    $self->{unigram} = LanguageModel::Table::PitmanYorProcess->new($self->{zerogram},
								   \&LanguageModel::Unigram::PitmanYor::reduceUnigram);
}

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->{unigram}->add($sentence->word($node), $NULL_CHAR);
	$self->addPosList($sentence->word($node), $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->{unigram}->remove($sentence->word($node), $NULL_CHAR);
	$self->removePosList($sentence->word($node), $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
    }
}

sub sampleSentence {
    my ($self, $sentence, $anneal) = @_;
    $self->LanguageModel::Unigram::TypeBased::sampleSentence($sentence, $anneal);
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    my ($leftCount, $rightCount, $segmentedList) = $self->getSegmentedList($region, $iter);
    my ($unsegmentedCount, $unsegmentedList) = $self->getUnsegmentedList($region, $iter);
    my $numState = scalar(@$segmentedList) + scalar(@$unsegmentedList);

    # use simpler sampling procedure if no other site is found
    return $self->sampleOneBoundary($region, $iter, $anneal) if ($numState <= 1);
    $self->{totalRuns}++;

    my ($proposal, $logQF, $logQB) = getFlipProposal(scalar(@$unsegmentedList), $numState);
    my ($splitList, $mergeList) = $self->setArrangements($proposal, $numState, $segmentedList, $unsegmentedList);

    my $llength = length($region->{left});
    my $currentLogProb = 0.0;
    foreach my $holder (@$splitList) {
	# remove unseg counts
	$self->{unigram}->remove($region->{unsegmented}, $NULL_CHAR);
	my $prob = $self->{unigram}->prob($region->{unsegmented}, $NULL_CHAR);
	$currentLogProb += log($prob);
    }
    foreach my $holder (@$mergeList) {
	# remove left and right counts
	$self->{unigram}->remove($region->{right}, $NULL_CHAR);
	my $prob = $self->{unigram}->prob($region->{right}, $NULL_CHAR);
	$self->{unigram}->remove($region->{left}, $NULL_CHAR);
	$prob *= $self->{unigram}->prob($region->{left}, $NULL_CHAR);
	$currentLogProb += log($prob);
    }

    my $proposalLogProb = 0.0;
    #  add unigram and unigram counts according to $splitList and $mergeList
    foreach my $holder (reverse(@$splitList)) {
	# add left and right counts
	my $prob = $self->{unigram}->prob($region->{left}, $NULL_CHAR);
	$self->{unigram}->add($region->{left}, $NULL_CHAR);
	$prob *= $self->{unigram}->prob($region->{right}, $NULL_CHAR);
	$self->{unigram}->add($region->{right}, $NULL_CHAR);
	$proposalLogProb += log($prob);
    }
    foreach my $holder (reverse(@$mergeList)) {
	# remove left and right counts
	my $prob = $self->{unigram}->prob($region->{unsegmented}, $NULL_CHAR);
	$self->{unigram}->add($region->{unsegmented}, $NULL_CHAR);
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
	    $self->{unigram}->remove($region->{right}, $NULL_CHAR);
	    $self->{unigram}->remove($region->{left}, $NULL_CHAR);
	    $self->{unigram}->add($region->{unsegmented}, $NULL_CHAR);
	}
	foreach my $holder (@$mergeList) {
	    # remove left and right counts
	    $self->{unigram}->remove($region->{unsegmented}, $NULL_CHAR);
	    $self->{unigram}->add($region->{left}, $NULL_CHAR);
	    $self->{unigram}->add($region->{right}, $NULL_CHAR);
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

    $self->{skipBlocks} = {};
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
