package LanguageModel::Bigram::Backward::TypeBased;

use strict;
use warnings;
use utf8;

use Carp::Assert;
use Scalar::Util qw/refaddr/;

use Text::Sentence;
use LanguageModel::Util qw/randPartitionLog getFlipProposal/;

use base qw/LanguageModel::Bigram::TypeBased/;

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->add($prev, $current); # reverse
	$self->addPosList($current, $sid, $iter->wpos, $sentence, $node) unless ($isFrozen);
	$prev = $current;
    }
    $self->{bigram}->add($prev, '$'); # reverse
}

sub removeSentence {
    my ($self, $sentence, $docID, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $sid = $sentence->{sid};
    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->remove($prev, $current); # reverse
	$self->removePosList($current, $sid, $iter->wpos) unless ($isFrozen);
	$prev = $current;
    }
    $self->{bigram}->remove($prev, '$'); # reverse
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
	$self->{bigram}->remove($region->{unsegmented}, $outRight);          # reverse
	my $prob = $self->{bigram}->prob($region->{unsegmented}, $outRight); # reverse
	$self->{bigram}->remove($outLeft, $region->{unsegmented});           # reverse
	$prob *= $self->{bigram}->prob($outLeft, $region->{unsegmented});    # reverse
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
	$self->{bigram}->remove($region->{right}, $outRight);              # reverse
	my $prob = $self->{bigram}->prob($region->{right}, $outRight);     # reverse
	$self->{bigram}->remove($region->{left}, $region->{right});        # reverse
	$prob *= $self->{bigram}->prob($region->{left}, $region->{right}); # reverse
	$self->{bigram}->remove($outLeft, $region->{left});                # reverse
	$prob *= $self->{bigram}->prob($outLeft, $region->{left});         # reverse
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
	my $prob = $self->{bigram}->prob($outLeft, $region->{left});       # reverse
	$self->{bigram}->add($outLeft, $region->{left});                   # reverse
	$prob *= $self->{bigram}->prob($region->{left}, $region->{right}); # reverse
	$self->{bigram}->add($region->{left}, $region->{right});           # reverse
	$prob *= $self->{bigram}->prob($region->{right}, $outRight);       # reverse
	$self->{bigram}->add($region->{right}, $outRight);                 # reverse
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
	my $prob = $self->{bigram}->prob($outLeft, $region->{unsegmented}); # reverse
	$self->{bigram}->add($outLeft, $region->{unsegmented});             # reverse
	$prob *= $self->{bigram}->prob($region->{unsegmented}, $outRight);  # reverse
	$self->{bigram}->add($region->{unsegmented}, $outRight);            # reverse
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
	    $self->{bigram}->remove($region->{right}, $outRight);       # reverse
	    $self->{bigram}->remove($region->{left}, $region->{right}); # reverse
	    $self->{bigram}->remove($outLeft, $region->{left});         # reverse
	    $self->{bigram}->add($outLeft, $region->{unsegmented});     # reverse
	    $self->{bigram}->add($region->{unsegmented}, $outRight);    # reverse
	}
	foreach my $holder (@$mergeList) {
	    # remove left and right counts
	    my $left = $holder->{node};
	    my $prevNode = $left->[Text::Sentence->PREV];
	    my $right = $left->[Text::Sentence->NEXT];
	    my $nextNode = $right->[Text::Sentence->NEXT];
	    my $outLeft = (defined($prevNode))? $holder->{sentence}->word($prevNode) : '$';
	    my $outRight = (defined($nextNode))? $holder->{sentence}->word($nextNode) : '$';
	    $self->{bigram}->remove($region->{unsegmented}, $outRight); # reverse
	    $self->{bigram}->remove($outLeft, $region->{unsegmented});  # reverse
	    $self->{bigram}->add($outLeft, $region->{left});            # reverse
	    $self->{bigram}->add($region->{left}, $region->{right});    # reverse
	    $self->{bigram}->add($region->{right}, $outRight);          # reverse
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
    my $newSegType = $self->LanguageModel::Bigram::Backward::sampleBoundary($region, $iter, $anneal);
    if ($newSegType == Text::Sentence->BOUNDARY) {
	$self->addPosList($region->{left}, $sid, $wpos, $sentence, $iter->{node});
	$self->addPosList($region->{right}, $sid, $wpos + $llength, $sentence, $iter->{node}->[Text::Sentence->NEXT]);
    } else {
	$self->addPosList($region->{unsegmented}, $sid, $wpos, $sentence, $iter->{node});
    }
    return $newSegType;
}

1;
