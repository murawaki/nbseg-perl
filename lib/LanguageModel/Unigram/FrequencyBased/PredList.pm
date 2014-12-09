package LanguageModel::Unigram::FrequencyBased::PredList;
#
# token-based sampler for text in which each sentence has a frequency count
#   incorporates a list of predicates
#
use strict;
use warnings;
use utf8;

use Carp::Assert;
use Math::Cephes qw/lgam psi/;
use List::Util qw/sum/;

use Text::Sentence;
use LanguageModel::Unigram::FrequencyBased;
use LanguageModel::Util qw/randPartitionLog logsumexpList/;

use base qw/LanguageModel::Unigram::FrequencyBased/;

# our $HP_DRAW_COUNT = 5;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	alpha => 0.1,
	a => 1, b => 1, # G(a,b), the prior for the hyperparameter
	zerogram => shift,
	opt => shift,
	unigram => {},
	predAlpha => 10,
	predList => {}, # pred -> { total: count, list: { head: count, ... } }
	tokenCount => 0,
    };
    $self->{alpha} = $self->{opt}->{alpha} if (defined($self->{opt}->{alpha}));
    $self->{opt}->{matchingSkip} = 0 unless (defined($self->{opt}->{matchingSkip}));
    $self->{opt}->{doSampleHyper} = 0 unless (defined($self->{opt}->{doSampleHyper}));
    bless($self, $class);
    return $self;
}

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;

    my $freq = $sentence->freq;
    my $iter = $sentence->nodeIterator;
    my $head;
    while ((my $node = $iter->next)) {
	$head = $sentence->word($node);
	$self->addWord($head, $freq);
    }
    while ((my ($pred, $freq) = each(%{$sentence->{predList}}))) {
	$self->addWord2Pred($pred, $head, $freq);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;

    my $freq = $sentence->freq;
    my $iter = $sentence->nodeIterator;
    my $head;
    while ((my $node = $iter->next)) {
	$head = $sentence->word($node);
	$self->removeWord($head, $freq);
    }
    while ((my ($pred, $freq) = each(%{$sentence->{predList}}))) {
	$self->removeWord2Pred($pred, $head, $freq);
    }
}

sub addWord2Pred {
    my ($self, $pred, $word, $freq) = @_;
    $freq = 1 unless (defined($freq));

    $self->{predList}->{$pred}->{total} += $freq;
    $self->{predList}->{$pred}->{list}->{$word} += $freq;
}

sub removeWord2Pred {
    my ($self, $pred, $word, $freq) = @_;
    $freq = 1 unless (defined($freq));

    my $predList = $self->{predList};
    $predList->{$pred}->{total} -= $freq;
    $predList->{$pred}->{list}->{$word} -= $freq;

    assert($predList->{$pred}->{total} >= 0);
    assert($predList->{$pred}->{list}->{$word} >= 0);
    if ($self->{predList}->{$pred}->{list}->{$word} <= 0) {
	delete($predList->{$pred}->{list}->{$word});
    }
    if ($predList->{$pred}->{total} <= 0) {
	assert(scalar(keys(%{$predList->{$pred}->{list}})) == 0);
	delete($predList->{$pred});
    }
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal, $freq) = @_;

    my $wpos = $iter->wpos;
    my $llength = length($region->{left});
    if ($region->{segType} == Text::Sentence->NONE_BOUNDARY) {
	$self->removeWord($region->{unsegmented}, $freq);
    } else {
	$self->removeWord($region->{left}, $freq);
	$self->removeWord($region->{right}, $freq);
    }
    if ($region->{isRightmost}) {
	my $head = ($region->{segType} == Text::Sentence->NONE_BOUNDARY)?
	    $region->{unsegmented} : $region->{right};
	my $sentence = $iter->{sentence};
	while ((my ($pred, $freq) = each(%{$sentence->{predList}}))) {
	    $self->removeWord2Pred($pred, $head, $freq);
	}
    }

    # frequency counts
    my $leftCount = $self->{unigram}->{$region->{left}} || 0;
    my $rightCount = $self->{unigram}->{$region->{right}} || 0;
    my $unsegmentedCount = $self->{unigram}->{$region->{unsegmented}} || 0;

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
	$logProbSeg += lgam($ALPHA * $self->zeroProb($region->{right}) + $rightCount + $freq);
	$logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{left}) + $leftCount);
	$logProbUnseg += lgam($ALPHA * $self->zeroProb($region->{right}) + $rightCount);
    }

    if ($region->{isRightmost}) {
	$logProbSeg += $self->predProb($region->{right}, $iter->{sentence}->{predList}); # , $self->zeroProb($region->{right}));
	$logProbUnseg += $self->predProb($region->{unsegmented}, $iter->{sentence}->{predList}); # , $self->zeroProb($region->{unsegmented}));
    }

    my $newNode = $iter->{node};
    my $newSegType = (randPartitionLog([$logProbSeg, $logProbUnseg], $anneal) == 0)?
	Text::Sentence->BOUNDARY : Text::Sentence->NONE_BOUNDARY;

    if ($region->{isRightmost}) {
	my $head = ($newSegType == Text::Sentence->NONE_BOUNDARY)?
	    $region->{unsegmented} : $region->{right};
	my $sentence = $iter->{sentence};
	while ((my ($pred, $freq) = each(%{$sentence->{predList}}))) {
	    $self->addWord2Pred($pred, $head, $freq);
	}
    }

    if ($newSegType == Text::Sentence->BOUNDARY) {
	if ($region->{segType} != Text::Sentence->BOUNDARY) {
	    $newNode = $iter->split($region->{left}, $region->{right});
	}
	$self->addWord($region->{left}, $freq);
	$self->addWord($region->{right}, $freq);
	return Text::Sentence->BOUNDARY;
    } else {
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->addWord($region->{unsegmented}, $freq);
	return Text::Sentence->NONE_BOUNDARY;
    }
}

sub predProb {
    my ($self, $head, $sentencePredList) = @_;
    # my ($self, $head, $sentencePredList, $backoff) = @_;

    # uniform prior over words
    my $predList = $self->{predList};
    # my $alpha = $self->{predAlpha};
    # my $alphaG0 = $alpha * $backoff;
    my $logProb = 0;
    while ((my ($pred, $freq) = each(%$sentencePredList))) {
	# ignore shared terms
	# $logProb += lgam((($predList->{$pred} && $predList->{$pred}->{list}->{$head}) || 0) + $freq + $alphaG0);
	$logProb += lgam((($predList->{$pred} && $predList->{$pred}->{list}->{$head}) || 0) + $freq);
    }
    return $logProb;
}


# sentence-based block sampling
# NOTE: we consider generation of a single sentence while we consider
# generation from multiple predicates
sub sampleSentenceBlock {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    $self->removeSentence($sentence);
    my $rawSentence = $sentence->asText(1);
    my $l = length($rawSentence);
    my $ALPHA = $self->{alpha};

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
	$constraints->{length($rawSentence)} = 1;
    }

    # forward-filtering
    my $logProbList = {}; # word -> prob
    my $alpha = [0]; # right position -> accumulated prob
    my $beta = [];   # right position -> word length - 1 -> accumulated prob
    for (my $t = 1; $t <= $l - 1; $t++) { # right position
	my $accumulatedList = [];
	for (my $k = $t - 1; $k >= 0; $k--) { # left position
	    my $word = substr($rawSentence, $k, $t - $k);
	    my $logProb = $logProbList->{$word} ||
		($logProbList->{$word} = log((($self->{unigram}->{$word} || 0) + $ALPHA * $self->zeroProb($word))
					     / ($self->{tokenCount} + $ALPHA)));
	    $beta->[$t]->[$t - $k - 1] = $logProb + $alpha->[$k];
	    push(@$accumulatedList, $logProb + $alpha->[$k]);
	    last if ($constraints->{$k}); # cannot go beyond a fixed boundary
	}
	$alpha->[$t] = (scalar(@$accumulatedList) > 1)? logsumexpList($accumulatedList) : $accumulatedList->[0];
    }
    {
	# for head, we condier generation from predicates
	my $accumulatedList = [];
	for (my $k = $l - 1; $k >= 0; $k--) { # left position
	    my $word = substr($rawSentence, $k, $l - $k);
	    my $logProb = $logProbList->{$word} ||
		($logProbList->{$word} = log((($self->{unigram}->{$word} || 0) + $ALPHA * $self->zeroProb($word))
					     / ($self->{tokenCount} + $ALPHA)));
	    $logProb += $self->predProb($word, $sentence->{predList}); # , $self->zeroProb($word));
	    $beta->[$l]->[$l - $k - 1] = $logProb + $alpha->[$k];
	    push(@$accumulatedList, $logProb + $alpha->[$k]);
	    last if ($constraints->{$k}); # cannot go beyond a fixed boundary
	}
    }

    # backward-sampling
    my $t = length($rawSentence);
    my $reverseWordList = [];
    while ($t > 0) {
	my $logProbList = [];
	my $indexList = [];
	for (my $i = 0; $i < scalar(@{$beta->[$t]}); $i ++) {
	    next unless (defined($beta->[$t]->[$i]));
	    push(@$logProbList, $beta->[$t]->[$i]);
	    push(@$indexList, $i);
	}
	my $idx = randPartitionLog($logProbList, $anneal);
	my $i = $indexList->[$idx];
	my $word = substr($rawSentence, $t - $i - 1, $i + 1);
	my $segType = ($constraints->{$t})? Text::Sentence->FIXED_BOUNDARY : Text::Sentence->BOUNDARY;
	push(@$reverseWordList, [$word, $segType]);
	$t -= $i + 1;
    }

    my @wordList = reverse(@$reverseWordList);
    $sentence->replaceFully(\@wordList);

    $self->addSentence($sentence);
}

sub sampleHyper {
    my ($self) = @_;
    $self->SUPER::sampleHyper;

    # # alpha is shared among predicates
    # my ($tableCount, $tokenCount) = (0, 0);
    # my $alpha = $self->{predAlpha};
    # while ((my ($pred, $struct) = each(%{$self->{predList}}))) {
    # 	while ((my ($head, $count) = each(%{$struct->{list}}))) {
    # 	    my $zeroProb = $alpha * $self->zeroProb($head);
    # 	    $tableCount += $zeroProb * (psi($zeroProb + $count) - psi($zeroProb));
    # 	    $tokenCount += $count;
    # 	}
    # }

    # # MCMC estimation
    # my $alphaOld = $self->{predAlpha};
    # my $alphaList = [];
    # for (1 .. $HP_DRAW_COUNT) {
    # 	$alphaOld = $self->drawAlphaFlat($alphaOld, $tableCount, $tokenCount);
    # 	push(@$alphaList, $alphaOld);
    # }
    # $self->{predAlpha} = sum(@$alphaList) / $HP_DRAW_COUNT;
}

sub calcLogProb {
    my ($self) = @_;

    # log likelihood of word generation
    my $ll = $self->SUPER::calcLogProb;

    # log likelihood of predicate-to-head generation
    # my $alpha = $self->{predAlpha};
    while ((my ($pred, $struct) = each(%{$self->{predList}}))) {
	$ll += lgam($alpha) - lgam($alpha + $struct->{total});
	while ((my ($head, $count) = each(%{$struct->{list}}))) {
	    # my $zeroProb = $alpha * $self->zeroProb($head);
	    # $ll += lgam($zeroProb + $count) - lgam($zeroProb);
	    $ll += lgam($count);
	}
    }
    return $ll;
}

1;
