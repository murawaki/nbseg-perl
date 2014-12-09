package LanguageModel::Bigram::Backward;

use strict;
use warnings;
use utf8;

use Math::GSL::Randist qw/gsl_ran_flat/;

use Text::Sentence;
use LanguageModel::Util qw/$rng randPartitionLog logsumexpList/;

use base qw/LanguageModel::Bigram/;

our $NULL_CHAR = '__NULL__';

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;

    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->add($prev, $current); # reverse
	$prev = $current;
    }
    $self->{bigram}->add($prev, '$'); # reverse
}

sub removeSentence {
    my ($self, $sentence, $docID, $isFrozen) = @_;

    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->remove($prev, $current); # reverse
	$prev = $current;
    }
    $self->{bigram}->remove($prev, '$'); # reverse
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    my $prevNode = $region->{node}->[Text::Sentence->PREV];
    my $outLeft = (defined($prevNode))? $iter->{sentence}->word($prevNode) : '$';
    my $nextNode = ($region->{segType} == Text::Sentence->NONE_BOUNDARY)?
	$region->{node}->[Text::Sentence->NEXT]
	: $region->{node}->[Text::Sentence->NEXT]->[Text::Sentence->NEXT];
    my $outRight = (defined($nextNode))? $iter->{sentence}->word($nextNode) : '$';
    if ($region->{segType} == Text::Sentence->NONE_BOUNDARY) {
	# NOTE: these operations decrement too much
	#   we can keep the unigram count of $outRight,
	#   but remove it for simplicity
	$self->{bigram}->remove($region->{unsegmented}, $outRight);  # reverse
	$self->{bigram}->remove($outLeft, $region->{unsegmented});   # reverse
    } else {
	$self->{bigram}->remove($region->{right}, $outRight);        # reverse
	$self->{bigram}->remove($region->{left}, $region->{right});  # reverse
	$self->{bigram}->remove($outLeft, $region->{left});          # reverse
    }

    # NOTE: approximation: do not fully update counts during sampling
    my $probSeg = $self->{bigram}->prob($region->{right}, $outRight) # reverse
	* $self->{bigram}->prob($region->{left}, $region->{right})   # reverse
	* $self->{bigram}->prob($outLeft, $region->{left});          # reverse
    my $probUnseg = $self->{bigram}->prob($region->{unsegmented}, $outRight) # revserse
	* $self->{bigram}->prob($outLeft, $region->{unsegmented});   # reverse

    if ($anneal != 1) {
	my $mass = $probSeg + $probUnseg;
	$probSeg = ($probSeg / $mass) ** $anneal;
	$probUnseg = ($probUnseg / $mass) ** $anneal;
    }

    my $newNode;
    if (gsl_ran_flat($rng->raw, 0, $probSeg + $probUnseg) < $probSeg) {
	if ($region->{segType} != Text::Sentence->BOUNDARY) {
	    $newNode = $iter->split($region->{left}, $region->{right});
	}
	$self->{bigram}->add($region->{right}, $outRight);       # reverse
	$self->{bigram}->add($region->{left}, $region->{right}); # reverse
	$self->{bigram}->add($outLeft, $region->{left});         # reverse
	return Text::Sentence->BOUNDARY;
    } else {
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->{bigram}->add($region->{unsegmented}, $outRight); # reverse
	$self->{bigram}->add($outLeft, $region->{unsegmented});  # reverse
	return Text::Sentence->NONE_BOUNDARY;
    }
}

# sentence-based block sampling
# see Mochihashi+ (2009)
sub sampleSentenceBlock {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    $self->removeSentence($sentence);
    my $rawSentence = $sentence->asText(1);
    my $L = length($rawSentence);

    my $constraints = {};
    {
	my $iter = $sentence->nodeIterator;
	my $wpos = 0;
	while ((my $node = $iter->next)) {
	    $wpos += length($sentence->word($node));
	    if ($sentence->segType($node) == Text::Sentence->FIXED_BOUNDARY) {
		$constraints->{$wpos} = 1;
	    }
	}
    }

    # forward-filtering
    my $logProbList = {}; # word -> next := prob
    my $alpha = []; $alpha->[$L]->[0] = 0; # right position -> word length - 1 := accumulated prob
    for (my $t = $L - 1;  $t >= 0; $t--) { # left position
	for (my $k = $t + 1; $k <= $L; $k++) { # right position
	    my $accumulatedList = [];
	    my $word = substr($rawSentence, $t, $k - $t);
	    for (my $j = $k + 1; $j <= $L + 1; $j++) { # right position of next word
		last if ($k < $L && $j > $L); # special case: '$'
		my $next = ($j <= $L)? substr($rawSentence, $k, $j - $k) : '$';
		my $logProb = ($logProbList->{$word}->{$next}
			       or ($logProbList->{$word}->{$next} = log($self->{bigram}->prob($word, $next))));
		push(@$accumulatedList, $logProb + $alpha->[$k]->[$j - $k - 1]);
		last if ($constraints->{$j}); # cannot go beyond a fixed boundary
	    }
	    $alpha->[$t]->[$k - $t -1] = (scalar(@$accumulatedList) > 1)? logsumexpList($accumulatedList) : $accumulatedList->[0];
	    last if ($constraints->{$k}); # cannot go beyond a fixed boundary
	}
    }

    # backward-sampling
    my $t = 0;
    my $word = '$';
    my $wordList = [];
    while ($t < $L) {
	my $i;
	if (scalar(@{$alpha->[$t]}) > 1) {
	    my $logMassList = [];
	    for (my $i = 0; $i < scalar(@{$alpha->[$t]}); $i++) { # length - 1
		my $next = substr($rawSentence, $t, $i + 1);
		my $logProb = ($logProbList->{$word}->{$next}
			       or ($logProbList->{$word}->{$next} = log($self->{bigram}->prob($word, $next))));
		push(@$logMassList, $logProb + $alpha->[$t]->[$i]);
	    }
	    $i = randPartitionLog($logMassList, $anneal);
	} else {
	    $i = 0;
	}
	$word = substr($rawSentence, $t, $i + 1);
	my $segType = ($constraints->{$t + $i + 1})? Text::Sentence->FIXED_BOUNDARY : Text::Sentence->BOUNDARY;
	push(@$wordList, [$word, $segType]);
	$t += $i + 1;
    }

    $sentence->replaceFully($wordList);
    $self->addSentence($sentence);
}

sub sampleLatent {
    my ($self, $sentence) = @_;

    my $iter = $sentence->nodeIterator;
    my $prev = '$';
    while ((my $node = $iter->next)) {
    	my $word = $sentence->word($node);
	$self->{bigram}->remove($prev, $word); # reverse
	$self->{bigram}->add($prev, $word);
	$prev = $word;
    }
    $self->{bigram}->remove('$', $prev);
    $self->{bigram}->add('$', $prev);
}

1;
