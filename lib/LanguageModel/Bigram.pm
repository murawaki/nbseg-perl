package LanguageModel::Bigram;
#
# bigram language model based on Hierarchical Dirichlet/Pitman-Yor processed
#
use strict;
use warnings;
use utf8;

use Math::GSL::Randist qw/gsl_ran_flat/;

use LanguageModel::Unigram;
use Text::Sentence;
use LanguageModel::Util qw/$rng randPartitionLog logsumexpList/;
use LanguageModel::Table::ZeroAdaptor;
use LanguageModel::Table::DirichletProcess; # explicit table tracking for hyperparameter estimation

our $NULL_CHAR = '__NULL__';

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	alpha0 => 10000, # 1000
	alpha1 => 10000, # 0.1,
	zerogramModel => shift,
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

    $self->{zerogram} = LanguageModel::Table::ZeroAdaptor->new($self->{zerogramModel}, !($self->{opt}->{nested}) );
    my $tableClass = ($self->{opt}->{type} eq 'Dirichlet')?
	'LanguageModel::Table::DirichletProcess' : 'LanguageModel::Table::PitmanYorProcess';
    $self->{unigram} = $tableClass->new($self->{zerogram}, \&reduceUnigram);
    $self->{bigram} = $tableClass->new($self->{unigram}, \&reduceBigram);
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

    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->add($current, $prev);
	$prev = $current;
    }
    $self->{bigram}->add('$', $prev);
}

sub removeSentence {
    my ($self, $sentence, $docID, $isFrozen) = @_;

    my $prev = '$';
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $current = $sentence->word($node);
	$self->{bigram}->remove($current, $prev);
	$prev = $current;
    }
    $self->{bigram}->remove('$', $prev);
}

sub calcLogProb {
    my ($self) = @_;

    return $self->{bigram}->logProb;
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
	$self->{bigram}->remove($region->{unsegmented}, $outLeft);
	$self->{bigram}->remove($outRight, $region->{unsegmented});
    } else {
	$self->{bigram}->remove($region->{left}, $outLeft);
	$self->{bigram}->remove($region->{right}, $region->{left});
	$self->{bigram}->remove($outRight, $region->{right});
    }

    # NOTE: approximation: do not update counts during sampling
    my $probSeg = $self->{bigram}->prob($region->{left}, $outLeft)
	* $self->{bigram}->prob($region->{right}, $region->{left})
	* $self->{bigram}->prob($outRight, $region->{right});
    my $probUnseg = $self->{bigram}->prob($region->{unsegmented}, $outLeft)
	* $self->{bigram}->prob($outRight, $region->{unsegmented});

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
	$self->{bigram}->add($region->{left}, $outLeft);
	$self->{bigram}->add($region->{right}, $region->{left});
	$self->{bigram}->add($outRight, $region->{right});
	return Text::Sentence->BOUNDARY;
    } else {
	if ($region->{segType} != Text::Sentence->NONE_BOUNDARY) {
	    $newNode = $iter->merge($region->{unsegmented});
	}
	$self->{bigram}->add($region->{unsegmented}, $outLeft);
	$self->{bigram}->add($outRight, $region->{unsegmented});
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
    my $logProbList = {}; # word -> prev := prob
    my $alpha = [[0]]; # right position -> word length - 1 := accumulated prob
    for (my $t = 1; $t <= $L; $t++) { # right position
	for (my $k = $t - 1; $k >= 0; $k--) { # left position
	    my $accumulatedList = [];
	    my $word = substr($rawSentence, $k, $t - $k);
	    for (my $j = $k - 1; $j >= -1; $j--) { # left position of previous word
		last if ($k > 0 && $j < 0); # special case: '$'
		my $prev = ($j >= 0)? substr($rawSentence, $j, $k - $j) : '$';
		my $logProb = ($logProbList->{$word}->{$prev}
			       or ($logProbList->{$word}->{$prev} = log($self->{bigram}->prob($word, $prev))));
		push(@$accumulatedList, $logProb + $alpha->[$k]->[$k - $j - 1]);
		last if ($constraints->{$j}); # cannot go beyond a fixed boundary
	    }
	    $alpha->[$t]->[$t - $k -1] = (scalar(@$accumulatedList) > 1)? logsumexpList($accumulatedList) : $accumulatedList->[0];
	    last if ($constraints->{$k}); # cannot go beyond a fixed boundary
	}
    }

    # backward-sampling
    my $t = $L;
    my $word = '$';
    my $reverseWordList = [];
    while ($t > 0) {
	my $i;
	if (scalar(@{$alpha->[$t]}) > 1) {
	    my $logMassList = [];
	    for (my $i = 0; $i < scalar(@{$alpha->[$t]}); $i++) { # length - 1
		my $prev = substr($rawSentence, $t - $i - 1, $i + 1);
		my $logProb = ($logProbList->{$word}->{$prev}
			       or ($logProbList->{$word}->{$prev} = log($self->{bigram}->prob($word, $prev))));
		push(@$logMassList, $logProb + $alpha->[$t]->[$i]);
	    }
	    $i = randPartitionLog($logMassList, $anneal);
	} else {
	    $i = 0;
	}
	$word = substr($rawSentence, $t - $i - 1, $i + 1);
	my $segType = ($constraints->{$t})? Text::Sentence->FIXED_BOUNDARY : Text::Sentence->BOUNDARY;
	push(@$reverseWordList, [$word, $segType]);
	$t -= $i + 1;
    }

    my @wordList = reverse(@$reverseWordList);
    $sentence->replaceFully(\@wordList);
    $self->addSentence($sentence);
}

sub sampleLatent {
    my ($self, $sentence) = @_;

    my $iter = $sentence->nodeIterator;
    my $prev = '$';
    while ((my $node = $iter->next)) {
    	my $word = $sentence->word($node);
	$self->{bigram}->remove($word, $prev);
	$self->{bigram}->add($word, $prev);
	$prev = $word;
    }
    $self->{bigram}->remove('$', $prev);
    $self->{bigram}->add('$', $prev);
}

sub sampleHyper {
    my ($self) = @_;
    # return unless($self->{opt}->{doSampleHyper});

    $self->{bigram}->sampleHyper;
    $self->{unigram}->sampleHyper;
}

sub setTitle {
    my ($self, $title) = @_;
    return $self->LanguageModel::Unigram::setTitle($title);
}

sub doSkip {
    my ($self, $region) = @_;
    $self->LanguageModel::Unigram::doSkip($region);
}

sub clearIteration {
    my ($self) = @_;
    $self->LanguageModel::Unigram::clearIteration;
}

sub reduceUnigram {
    my ($current, $context) = @_;
    return ($current, $NULL_CHAR);
}

sub reduceBigram {
    my ($current, $context) = @_;
    return ($current, $NULL_CHAR);
}

1;
