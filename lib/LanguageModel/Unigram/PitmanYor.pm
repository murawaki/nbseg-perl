package LanguageModel::Unigram::PitmanYor;
#
# the unigram language model based on a Pitman-Yor process
#   with sentence-based block sampling (no token-based sampling)
#
use strict;
use warnings;
use utf8;

use Text::Sentence;
use LanguageModel::Unigram;
use LanguageModel::Util qw/$rng randPartitionLog logsumexpList/;
use LanguageModel::Table::PitmanYorProcess;
use Math::GSL::Randist qw/gsl_ran_flat/;

use base qw/LanguageModel::Unigram/;

our $NULL_CHAR = '__NULL__';

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	zerogramModel => shift,
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

    # nested -> not static
    $self->{zerogram} = LanguageModel::Table::ZeroAdaptor->new($self->{zerogramModel}, !($self->{opt}->{nested}) );
    $self->{unigram} = LanguageModel::Table::PitmanYorProcess->new($self->{zerogram}, \&reduceUnigram);
}

sub addSentence {
    my ($self, $sentence, $isFrozen) = @_;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->{unigram}->add($sentence->word($node), $NULL_CHAR);
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->{unigram}->remove($sentence->word($node), $NULL_CHAR);
    }
}

# sub sampleSentence {
#     die("token-based sampling is not supported\n");
# }

sub sampleBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    if ($region->{segType} == Text::Sentence->NONE_BOUNDARY) {
	$self->{unigram}->remove($region->{unsegmented}, $NULL_CHAR);
    } else {
	$self->{unigram}->remove($region->{left}, $NULL_CHAR);
	$self->{unigram}->remove($region->{right}, $NULL_CHAR);
    }
    my $probUnseg = $self->{unigram}->prob($region->{unsegmented}, $NULL_CHAR);
    my $probSeg = $self->{unigram}->prob($region->{left}, $NULL_CHAR)
	* $self->{unigram}->prob($region->{right}, $NULL_CHAR);
    if ($anneal != 1) {
	my $mass = $probSeg + $probUnseg;
	$probSeg = ($probSeg / $mass) ** $anneal;
	$probUnseg = ($probUnseg / $mass) ** $anneal;
    }

    my $newSegType;
    my $newNode = $iter->{node};
    if (gsl_ran_flat($rng->raw, 0, $probSeg + $probUnseg) < $probSeg) {
	if ($region->{segType} != Text::Sentence->BOUNDARY) {
	    $newNode = $iter->split($region->{left}, $region->{right});
	}
	$self->{unigram}->add($region->{left}, $NULL_CHAR);
	$self->{unigram}->add($region->{right}, $NULL_CHAR);
	return Text::Sentence->BOUNDARY;
    } else {
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->{unigram}->add($region->{unsegmented}, $NULL_CHAR);
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
    }

    # forward-filtering
    my $logProbList = {}; # word := prob
    my $alpha = [0]; # right position := accumulated prob
    my $beta = [];   # right position -> word length - 1 := accumulated prob
    for (my $t = 1; $t <= $L; $t++) { # right position
	my $accumulatedList = [];
	for (my $k = $t - 1; $k >= 0; $k--) { # left position
	    my $word = substr($rawSentence, $k, $t - $k);
	    my $logProb = $logProbList->{$word} ||
		($logProbList->{$word} = log($self->{unigram}->prob($word, $NULL_CHAR)));
	    $beta->[$t]->[$t - $k - 1] = $logProb + $alpha->[$k];
	    push(@$accumulatedList, $logProb + $alpha->[$k]);
	    last if ($constraints->{$k}); # cannot go beyond a fixed boundary
	}
	$alpha->[$t] = (scalar(@$accumulatedList) > 1)? logsumexpList($accumulatedList) : $accumulatedList->[0];
    }

    # backward-sampling
    my $t = $L;
    my $reverseWordList = [];
    while ($t > 0) {
	my $i = (scalar(@{$beta->[$t]}) > 1)? randPartitionLog($beta->[$t], $anneal) : 0;
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
    $self->{unigram}->sampleHyper;
}

sub calcLogProb {
    my ($self) = @_;
    return $self->{unigram}->logProb;
}

sub reduceUnigram {
    my ($current, $context) = @_;
    return ($current, $NULL_CHAR);
}

1;
