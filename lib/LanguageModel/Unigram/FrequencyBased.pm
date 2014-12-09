package LanguageModel::Unigram::FrequencyBased;
#
# token-based sampler for text in which each sentence has a frequency count
#
use strict;
use warnings;
use utf8;

use Carp::Assert;
use Math::Cephes qw/lgam/;

use Text::Sentence;
use LanguageModel::Unigram;
use LanguageModel::Util qw/randPartitionLog/;

use base qw/LanguageModel::Unigram/;

our $HP_DRAW_COUNT = 5;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	alpha => 0.1,
	a => 1, b => 1, # G(a,b), the prior for the hyperparameter
	zerogram => shift,
	opt => shift,
	unigram => {},
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
    while ((my $node = $iter->next)) {
	$self->addWord($sentence->word($node), $freq);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;

    my $freq = $sentence->freq;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->removeWord($sentence->word($node), $freq);
    }
}

sub sampleSentence {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    my $freq = $sentence->freq;
    my $iter = $sentence->charIterator;
    while ((my $region = $iter->next)) {
	next if ($self->doSkip($region));
	my $segType = $self->sampleBoundary($region, $iter, $anneal, $freq);
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

    my $newNode = $iter->{node};
    if (randPartitionLog([$logProbSeg, $logProbUnseg], $anneal) == 0) {
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

1;
