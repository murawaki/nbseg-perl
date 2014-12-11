package LanguageModel::Zerogram;
#
# (katakana) zerogram or character n-gram model
#   using hierarchical Dirichlet/Pitman-Yor process
#
use strict;
use warnings;
use utf8;

use JSON;
use Math::GSL::RNG qw/gsl_rng_uniform_int/;

use LanguageModel::Util qw/$rng/;
use LanguageModel::Table::DirichletProcess;
use LanguageModel::Table::PitmanYorProcess;
use LanguageModel::Table::Root;

our $INV_SIZE_OF_CHARS = 1.0 / 88; # 0x30A1 - 0x30F6 + 0x30FC + '#'
our $NULL_CHAR = '__NULL__';
our $POOL_LIMIT = 100000;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	opt => shift,
	alpha0 => 60.0, alpha1 => 15.0,
	d1 => 0.3, theta1 => 4.0, d0 => 0.1, theta0 => 40.0,
	pool => {},
	poolLimit => $POOL_LIMIT,
	keyCount => 0,
	time => 0,
    };
    $self->{opt}->{type} = 'Dirichlet' unless (defined($self->{opt}->{type}));
    $self->{opt}->{doSampleHyper} = 0 unless (defined($self->{opt}->{doSampleHyper}));
    bless($self, $class);

    my $invSizeOfChars = $self->{opt}->{invSizeOfChars} || $INV_SIZE_OF_CHARS;
    $self->{zerogram} = LanguageModel::Table::Root->new($invSizeOfChars, \&drawChar);
    my $tableClass = ($self->{opt}->{type} eq 'Dirichlet')?
	'LanguageModel::Table::DirichletProcess' : 'LanguageModel::Table::PitmanYorProcess';
    $self->{unigram} = $tableClass->new($self->{zerogram}, \&reduceUnigram);
    $self->{bigram} = $tableClass->new($self->{unigram}, \&reduceBigram);
    if ($self->{opt}->{type} eq 'Dirichlet') {
	$self->{unigram}->alpha($self->{alpha0});
	$self->{bigram}->alpha($self->{alpha1});
    } else {
	$self->{unigram}->theta($self->{theta0});
	$self->{bigram}->theta($self->{theta1});
	$self->{unigram}->d($self->{d0});
	$self->{bigram}->d($self->{d1});
    }
    return $self;
}

sub addWord {
    my ($self, $name) = @_;

    my $p = '#';
    foreach my $c (split(//, $name)) {
	$self->{bigram}->add($c, $p);
	$p = $c;
    }
    $self->{bigram}->add('#', $p);
    $self->{pool} = {};
    $self->{keyCount} = 0;
}

sub removeWord {
    my ($self, $name) = @_;

    my $p = '#';
    foreach my $c (split(//, $name)) {
	$self->{bigram}->remove($c, $p);
	$p = $c;
    }
    $self->{bigram}->remove('#', $p);
    $self->{pool} = {};
    $self->{keyCount} = 0;
}

sub drawWord {
    my ($self) = @_;

    my $word = '';
    my $p = '#';
    my $bigram = $self->{bigram};
    while (1) {
	my $c = $bigram->draw($p);
	return $word if ($c eq '#');
	$word .= $c;
	$p = $c;
    }
}

sub getProb {
    my ($self, $name) = @_;

    $self->{time}++;
    my $record = $self->{pool}->{$name};
    if (defined ($record)) {
	$record->[1] = $self->{time};
	return $record->[0];
    } else {
	my $prob = $self->calcProb($name);
	$self->{pool}->{$name} = [$prob, $self->{time}];
	if (++$self->{keyCount} >= $self->{poolLimit}) {
	    $self->gc;
	}
	return $prob;
    }
}

sub gc {
    my ($self) = @_;

    my $pool = $self->{pool};
    my @sorted = sort { $pool->{$a}->[1] <=> $pool->{$b}->[1] } keys(%$pool);

    # printf STDERR ("running GC\tbefore: %d\n", scalar(@sorted));
    my $newLimit = int($self->{poolLimit} * 0.3);
    for (my $i = 0; $i < $newLimit; $i++) {
	delete($pool->{$sorted[$i]});
    }
    $self->{keyCount} = scalar(keys(%$pool));
    # printf STDERR ("running GC\tafter: %d\n", $self->{keyCount});
}

sub calcProb {
    my ($self, $name) = @_;

    my $prob = 1.0;
    my $p = '#';
    my $bigram = $self->{bigram};
    foreach my $c (split(//, $name)) {
	$prob *= $bigram->prob($c, $p);
	$p = $c;
    }
    $prob *= $bigram->prob('#', $p);
    return ($prob > 1e-95)? $prob : 1e-95; # check underflow
}

sub calcLogProb {
    my ($self) = @_;
    return $self->{bigram}->logProb;
}

sub sampleHyper {
    my ($self) = @_;

    # return unless($self->{opt}->{doSampleHyper});
    $self->{bigram}->sampleHyper;
    $self->{unigram}->sampleHyper;
}

sub dump {
    my ($self, $stream) = @_;

    my $json = JSON->new->allow_nonref->allow_blessed->convert_blessed;
    $stream->print($json->encode($self));
}

sub TO_JSON {
    my ($self) = @_;

    # remove recursive links, refs to subroutines
    return {
	opt => $self->{opt},
	zerogram => $self->{zerogram},
	unigram => $self->{unigram},
	bigram => $self->{bigram},
    };
}

sub load {
    my ($this, $stream) = @_;
    my $class = ref($this) || $this;

    my $buf = '';
    while ((my $line = $stream->getline)) {
	$buf .= $line;
    }
    my $json = JSON->new->allow_nonref;
    my $self = $json->decode($buf);
    bless($self, $class);
    $self->{pool} = {};
    $self->{poolLimit} = $POOL_LIMIT;
    $self->{keyCount} = 0;
    $self->{time} = 0;

    my $tableClass = 'LanguageModel::Table::' .
	(($self->{opt}->{type} eq 'Dirichlet')? 'DirichletProcess' : 'PitmanYorProcess');
    $tableClass->load($self->{unigram}, $self->{zerogram}, \&reduceUnigram);
    $tableClass->load($self->{bigram}, $self->{unigram}, \&reduceBigram);
    LanguageModel::Table::Root->load($self->{zerogram}, \&drawChar);
    return $self;
}

sub reduceUnigram {
    my ($current, $context) = @_;
    return ($NULL_CHAR, $NULL_CHAR);
}

sub reduceBigram {
    my ($current, $context) = @_;
    return ($current, $NULL_CHAR);
}

sub drawChar {
    # my ($context) = @_;
    my $r = gsl_rng_uniform_int($rng->raw, 87); # 0x30A1 - 0x30F6 + 0x30FC - + '#' - 1
    if ($r >= 87) {
	return '#';
    } elsif ($r >= 86) {
	return 'ー'; # 0x30FC
    } else {
	return chr(0x30A1 + $r);
    }
}


1;
