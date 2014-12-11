#!/bin/env perl
#
# segment sentences with the language model
# input/output format: "sentence" for each line
#
use strict;
use warnings;
use utf8;

use Getopt::Long;
use IO::File;

use LanguageModel::Zerogram;
use LanguageModel::Zerogram::UniChar;
use LanguageModel::Unigram;
use LanguageModel::Unigram::PitmanYor;
use LanguageModel::Unigram::TypeBased;
use LanguageModel::Unigram::PitmanYor::TypeBased;
use LanguageModel::Bigram;
use LanguageModel::Bigram::TypeBased;
use LanguageModel::Bigram::Predicate;
use LanguageModel::Bigram::Predicate::TypeBased;
use LanguageModel::Bigram::Backward;
use LanguageModel::Bigram::Backward::TypeBased;
use LanguageModel::Bigram::Backward::Predicate;
use LanguageModel::Bigram::Backward::Predicate::TypeBased;
use LanguageModel::Util qw/shuffle/;
use Text::Sentence;

binmode(STDIN,  ':utf8');
binmode(STDOUT, ':utf8');
binmode(STDERR, ':utf8');

my $opt = { type => 'Dirichlet',
	    iter => 20, randInit => 0.8, typeLevel => 0, sampleHyper => 0,
	    annealStart => 1, annealStop => 1, annealIter => 0,
	    nested => 0,
	    sentenceBlockIter => 0 };

sub main {
    GetOptions($opt,
	       'input=s',       # list of raw sentences
	       'type=s',        # Dirichlet or PitmanYor
	       'zero=s',        # zerogram model
	       'unichar=f',     # unichar zerogram model (specify boundary prob)
	       'charsize=i',    # number of characters for zerogram
	       'reference=s',   # goldstandard segmentation
	       'iter=i',        # number of iterations
	       'randInit=f',    # initialize segmentation by putting boundaries with the specified probability
	       'seed=i',        # random seed (default: time)
	       'bigram',        # use bigram model instad of unigram
	       'backward',      # backward transition (bigram model only)
	       'predicate',     # use predicate instead of EOS (bigram model only)
	       'typeLevel',     # type-based block sampling
	       'sentenceBlockIter=i', # sentence-based block sampling for the first N iterations
	       'sampleHyper',   # perform hyperparameter sampling
	       'nested',        # treat a zerogram model as a nested HPY
	                        # NOTE: specifying both 'reference' and 'zero' cause double counting
	       'saveZero=s',    # save the zerogram model to a specified file
	       'debug',         # pring debug messages
	       'annealStart=f', # start at this annealing temperature (usually less than 1.0)
	       'annealStop=f',  # stop at this annealing temperature (very large for MAP estimation)
	       'annealIter=i',  # do annealing for the first specified iterations
	);
    &LanguageModel::Util::setSeed($opt->{seed}) if (defined($opt->{seed}));
    my $zero;
    if ($opt->{zero} || $opt->{nested}) {
	$zero = &loadZerogram($opt->{zero}, $opt->{nested});
    } else {
	$zero = &loadUniChar($opt->{input}, $opt->{unichar});
    }

    my $modelClass = 'LanguageModel::' . (($opt->{bigram})? 'Bigram' : 'Unigram');
    my $modelOpt = { nested => $opt->{nested} };
    if ($opt->{type} eq 'PitmanYor') {
	if ($opt->{bigram}) {
	    $modelOpt->{type} = 'PitmanYor';
	} else {
	    $modelClass .= '::PitmanYor';
	}
    }
    $modelClass .= ($opt->{backward})? '::Backward' : '';
    $modelClass .= ($opt->{predicate})? '::Predicate' : '';
    $modelClass .= ($opt->{typeLevel})? '::TypeBased' : '';
    printf STDERR ("model %s (%s)\n", $modelClass, $opt->{type}) if ($opt->{debug});
    my $model = $modelClass->new($zero, $modelOpt);

    &loadInput($opt->{reference}, $model, 1) if ($opt->{reference});
    my $compoundList = &loadInput($opt->{input}, $model);
    printf STDERR ("logprob\t%f\n", $model->calcLogProb) if ($opt->{debug});

    my $annealLevel = 1;
    for my $iter (1 .. $opt->{iter}) {
	printf STDERR ("iteration %d\n", $iter) if ($opt->{debug});
	if ($opt->{annealIter} > 0) {
	    if ($iter < $opt->{annealIter}) {
		$annealLevel = $opt->{annealStart} * (($opt->{annealStop} / $opt->{annealStart}) ** (($iter - 1) / $opt->{annealIter}));
	    } else {
		$annealLevel = $opt->{annealStop};
	    }
	    printf STDERR ("annealing level: %f\n", $annealLevel);
	} else {
	    $annealLevel = 1;
	}

	shuffle($compoundList);
	if ($opt->{sentenceBlockIter}-- > 0) {
	    foreach my $compound (@$compoundList) {
		$model->sampleSentenceBlock($compound, $annealLevel);
	    }
	} else {
	    foreach my $compound (@$compoundList) {
		$model->sampleSentence($compound, $annealLevel);
	    }
	}
	&sampleHyper($model, $zero) if ($opt->{sampleHyper});
	printf STDERR ("logprob\t%f\n", $model->calcLogProb) if ($opt->{debug});
	$model->clearIteration;
    }

    foreach my $compound (sort { $a->{sid} <=> $b->{sid} } (@$compoundList)) {
	printf("%s\n", $compound->asText);
    }

    if ($opt->{saveZero}) {
	printf STDERR ("save model to %s...\n", $opt->{saveZero}) if ($opt->{debug});
	my $f = IO::File->new($opt->{saveZero}, 'w') or die;
	$f->binmode(':utf8');
	$model->{zerogramModel}->dump($f);
	$f->close;
    }
}

sub loadZerogram {
    my ($path, $isNested) = @_;

    if (defined($path)) {
	my $f = IO::File->new($path) or die;
	$f->binmode(':utf8');
	my $zero = LanguageModel::Zerogram->load($f);
	$f->close;
	return $zero;
    } else {
	die("load existing zerogram or use a nested model") unless ($isNested);
	my $zeroOpt = { type => $opt->{type} };
	if ($opt->{charsize}) {
	    $zeroOpt->{invSizeOfChars} = 1.0 / $opt->{charsize};
	    printf STDERR ("inv char size: %f\n", $zeroOpt->{invSizeOfChars});
	}
	return LanguageModel::Zerogram->new($zeroOpt);
    }
}

sub loadUniChar {
    my ($path, $ps) = @_;

    my $zero = LanguageModel::Zerogram::UniChar->new({ ps => $ps });
    my $f = IO::File->new($path) or die;
    $f->binmode(':utf8');
    while ((my $input = $f->getline)) {
	chomp($input);
	my @list = split(/\t/, $input);
	my $textWord = shift(@list);
	$zero->addWord($textWord);
    }
    $f->close;
    return $zero;
}

sub loadInput {
    my ($path, $model, $isFrozen) = @_;
    $isFrozen = 0 unless (defined($isFrozen));

    my $count = 0;
    my $compoundList = [];
    my $f = IO::File->new($path) or die;
    $f->binmode(':utf8');
  outer:
    while ((my $input = $f->getline)) {
	chomp($input);
	my @list = split(/\t/, $input);
	my $textWord = shift(@list);
	$textWord = '+' . $textWord . '+' unless ($textWord =~ /^\+.*\+$/);
	my $compound = Text::Sentence->createFromSegmentedText($textWord);
	$compound->randomize($opt->{randInit}) if (!$isFrozen && $opt->{randInit} >= 0);
	$compound->{sid} = $count++;
	if ($opt->{predicate}) {
	    # my $predList = {};
	    while ((my $pred = shift(@list))) {
		next if ($pred =~ /\:連体$/); # HACK!
		# $predList->{$pred} += 1;
		$compound->{pred} = $pred;
	    }
	    # $compound->{predList} = $predList;
	}
	$model->addSentence($compound, $isFrozen);
	push(@$compoundList, $compound);
    }
    return $compoundList;
}

sub sampleHyper {
    my ($model, $zero) = @_;

    if ($opt->{type} eq 'Dirichlet') {
	if ($opt->{bigram}) {
	    my $bAlphaOld = $model->{bigram}->alpha;
	    my $uAlphaOld = $model->{unigram}->alpha;
	    $model->sampleHyper;
	    my $bAlphaNew = $model->{bigram}->alpha;
	    my $uAlphaNew = $model->{unigram}->alpha;
	    printf STDERR ("bigram alpha %f -> %f\n", $bAlphaOld, $bAlphaNew) if ($opt->{debug});
	    printf STDERR ("unigram alpha %f -> %f\n", $uAlphaOld, $uAlphaNew) if ($opt->{debug});
	} else {
	    my $alphaOld = $model->{alpha};
	    my $predAlphaOld = $model->{predAlpha};
	    $model->sampleHyper;
	    my $alphaNew = $model->{alpha};
	    my $predAlphaNew = $model->{predAlpha};
	    printf STDERR ("alpha %f -> %f\n", $alphaOld, $alphaNew) if ($opt->{debug});
	    printf STDERR ("pred alpha %f -> %f\n", $predAlphaOld, $predAlphaNew) if ($opt->{debug} && $opt->{pred});
	}
    } else { # PitmanYor
	if ($opt->{bigram}) {
	    my $b_dOld = $model->{bigram}->{d};
	    my $b_thetaOld = $model->{bigram}->{theta};
	    my $u_dOld = $model->{unigram}->{d};
	    my $u_thetaOld = $model->{unigram}->{theta};
	    $model->sampleHyper;
	    my $b_dNew = $model->{bigram}->{d};
	    my $b_thetaNew = $model->{bigram}->{theta};
	    my $u_dNew = $model->{unigram}->{d};
	    my $u_thetaNew = $model->{unigram}->{theta};
	    if ($opt->{debug}) {
		printf STDERR ("bigram d %f -> %f\n", $b_dOld, $b_dNew);
		printf STDERR ("bigram theta %f -> %f\n", $b_thetaOld, $b_thetaNew);
		printf STDERR ("unigram d %f -> %f\n", $u_dOld, $u_dNew);
		printf STDERR ("unigram theta %f -> %f\n", $u_thetaOld, $u_thetaNew);
	    }
	} else {
	    my $dOld = $model->{unigram}->{d};
	    my $thetaOld = $model->{unigram}->{theta};
	    $model->sampleHyper;
	    my $dNew = $model->{unigram}->{d};
	    my $thetaNew = $model->{unigram}->{theta};
	    printf STDERR ("d %f -> %f\n", $dOld, $dNew) if ($opt->{debug});
	    printf STDERR ("theta %f -> %f\n", $thetaOld, $thetaNew) if ($opt->{debug});
	}
    }
    if ($opt->{nested}) {
	if ($opt->{type} eq 'Dirichlet') {
	    my $b_alphaOld = $zero->{bigram}->{alpha};
	    my $u_alphaOld = $zero->{unigram}->{alpha};
	    $zero->sampleHyper;
	    my $b_alphaNew = $zero->{bigram}->{alpha};
	    my $u_alphaNew = $zero->{unigram}->{alpha};
	    if ($opt->{debug}) {
		printf STDERR ("char bigram alpha %f -> %f\n", $b_alphaOld, $b_alphaNew);
		printf STDERR ("char unigram alpha %f -> %f\n", $u_alphaOld, $u_alphaNew);
	    }
	} else {
	    my $b_dOld = $zero->{bigram}->{d};
	    my $b_thetaOld = $zero->{bigram}->{theta};
	    my $u_dOld = $zero->{unigram}->{d};
	    my $u_thetaOld = $zero->{unigram}->{theta};
	    $zero->sampleHyper;
	    my $b_dNew = $zero->{bigram}->{d};
	    my $b_thetaNew = $zero->{bigram}->{theta};
	    my $u_dNew = $zero->{unigram}->{d};
	    my $u_thetaNew = $zero->{unigram}->{theta};
	    if ($opt->{debug}) {
		printf STDERR ("char bigram d %f -> %f\n", $b_dOld, $b_dNew);
		printf STDERR ("char bigram theta %f -> %f\n", $b_thetaOld, $b_thetaNew);
		printf STDERR ("char unigram d %f -> %f\n", $u_dOld, $u_dNew);
		printf STDERR ("char unigram theta %f -> %f\n", $u_thetaOld, $u_thetaNew);
	    }
	}
    }
}

unless (caller) {
    &main;
}

1;
