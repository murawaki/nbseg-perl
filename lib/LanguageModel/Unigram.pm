package LanguageModel::Unigram;
#
# Dirichlet process unigram language model
#   with token-based sampling and sentence-based block sampling
# this does not use LanguageModel::Table because table tracking is uncessesary
#
use strict;
use warnings;
use utf8;

use Carp::Assert;
use List::Util qw/max sum/;
use Regexp::Assemble;
use Math::Cephes qw/lgam psi/;
use Math::GSL::Randist qw/gsl_ran_flat/;

use Text::Sentence;
use LanguageModel::Util qw/$rng randPartitionLog logsumexpList/;
use LanguageModel::Table::DirichletProcess;

use base qw/LanguageModel::Table::DirichletProcess/; # for hyperparameter estimation

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
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->addWord($sentence->word($node));
    }
}

sub removeSentence {
    my ($self, $sentence, $isFrozen) = @_;
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	$self->removeWord($sentence->word($node));
    }
}

sub addWord {
    my ($self, $name, $freq) = @_;
    $freq = 1 unless(defined($freq));

    $self->{tokenCount} += $freq;;
    {
	no warnings qw/uninitialized/;
	return $self->{unigram}->{$name} += $freq;
    };
}

sub removeWord {
    my ($self, $name, $freq) = @_;
    $freq = 1 unless(defined($freq));

    $self->{tokenCount} -= $freq;;
    my $count = ($self->{unigram}->{$name} -= $freq);
    if ($count <= 0) {
	delete($self->{unigram}->{$name});
	assert($count >= 0, sprintf("no such word found: %s", $name)) if (DEBUG);
    }
    return $count;
}

sub zeroProb {
    my ($self, $name) = @_;
    return $self->{zerogram}->getProb($name);
}

sub sampleSentence {
    my ($self, $sentence, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    my $iter = $sentence->charIterator;
    while ((my $region = $iter->next)) {
	next if ($self->doSkip($region));
	my $segType = $self->sampleBoundary($region, $iter, $anneal);
    }
}

sub sampleBoundary {
    my ($self, $region, $iter, $anneal) = @_;

    if ($region->{segType} == Text::Sentence->NONE_BOUNDARY) {
	$self->removeWord($region->{unsegmented});
    } else {
	$self->removeWord($region->{left});
	$self->removeWord($region->{right});
    }
    # frequency counts
    my $left = $self->{unigram}->{$region->{left}} || 0;
    my $right = $self->{unigram}->{$region->{right}} || 0;
    my $unsegmented = $self->{unigram}->{$region->{unsegmented}} || 0;

    my $ALPHA = $self->{alpha};
    my $probUnseg = ($unsegmented + $ALPHA * $self->zeroProb($region->{unsegmented}));
	# / ($self->{tokenCount} + $ALPHA); # ignore shared term
    my $probSeg = ($left + $ALPHA * $self->zeroProb($region->{left}));
	# / ($self->{tokenCount} + $ALPHA); # ignore shared term
    $probSeg *= ($right + ($region->{left} eq $region->{right}? 1 : 0) + $ALPHA * $self->zeroProb($region->{right}))
	/ ($self->{tokenCount} + 1 + $ALPHA);
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
	$self->addWord($region->{left});
	$self->addWord($region->{right});
	return Text::Sentence->BOUNDARY;
    } else {
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->addWord($region->{unsegmented});
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
		($logProbList->{$word} = log((($self->{unigram}->{$word} || 0) + $ALPHA * $self->zeroProb($word))
					     / ($self->{tokenCount} + $ALPHA)));
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

    # NOTE: Metropolis-Hastings correction step is necessary here in
    #   order to correctly sample from the posterior distribution

    my @wordList = reverse(@$reverseWordList);
    $sentence->replaceFully(\@wordList);
    $self->addSentence($sentence);
}

sub sampleHyper {
    my ($self) = @_;

    # sampling of concentration parameter alpha
    # based on (Escobar and West, 1995)
    # we need:
    # * table count: k
    # * token count: n
    #
    # we first draw eta from a beta distribution
    # we then draw alpha from a mixture of gammas
    my ($k, $n) = $self->estimateTableCount;

    # MCMC estimation
    my $alphaOld = $self->{alpha};
    my $alphaList = [];
    for (1 .. $HP_DRAW_COUNT) {
	$alphaOld = $self->drawAlphaFlat($alphaOld, $k, $n);
	push(@$alphaList, $alphaOld);
    }
    $self->{alpha} = sum(@$alphaList) / $HP_DRAW_COUNT;
}

sub estimateTableCount {
    my ($self) = @_;

    # NOTE: this is heavy
    # approximating table counts
    # See Blunsom et al. 2009
    # E[t_w] = alpha * P(w) * (digamma(P(w) + n_w)) −  (digamma(P(w)))
    my ($tableCount, $tokenCount) = (0, 0);
    my $alpha = $self->{alpha};
    while ((my ($name, $count) = each(%{$self->{unigram}}))) {
	my $zeroProb = $alpha * $self->zeroProb($name);
	$tableCount += $zeroProb * (psi($zeroProb + $count) - psi($zeroProb));
	$tokenCount += $count;
    }
    return wantarray? ($tableCount, $tokenCount) : $tableCount;
}

sub calcLogProb {
    my ($self) = @_;

    my $alpha = $self->{alpha};
    my $ll = lgam($alpha) - lgam($alpha + $self->{tokenCount});
    while ((my ($name, $count) = each(%{$self->{unigram}}))) {
	my $zeroProb = $alpha * $self->zeroProb($name);
	$ll += lgam($zeroProb + $count) - lgam($zeroProb);
    }
    return $ll;
}

# Matching skip:
#   heuristic skip approximation
#
#   perform sampling only when the corresponding region contains a
#   substring of the noun phrase
#
sub setTitle {
    my ($self, $title) = @_;
    return unless ($self->{opt}->{matchingSkip});

    my $ra = Regexp::Assemble->new;
    foreach my $c (split(//, $title)) {
	$ra->add($c);
    }
    $self->{titleRE} = $ra->re;
    $self->{positionCount} = $self->{skipCount} = 0;
}

sub doSkip {
    my ($self, $region) = @_;
    return 0 unless ($self->{opt}->{matchingSkip});

    $self->{positionCount}++;
    if ($region->{unsegmented} =~ /$self->{titleRE}/) {
	return 0;
    } else {
	$self->{skipCount}++;
	return 1;
    }
}

sub clearIteration {
    my ($self) = @_;
    $self->{positionCount} = $self->{skipCount} = 0;
}

1;
