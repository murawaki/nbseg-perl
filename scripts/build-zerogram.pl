#!/bin/env perl
#
# input file format: [word] [freq] for each line
#
use strict;
use warnings;
use utf8;

use Getopt::Long;
use IO::File;
use LanguageModel::Util qw/shuffle/; # avoid List::Util::shuffle as it creates a new array
use LanguageModel::Zerogram;

binmode(STDIN,  ':utf8');
binmode(STDOUT, ':utf8');
binmode(STDERR, ':utf8');

my $opt = { type => 'Dirichlet', iter => 2, debug => 0 };

sub main {
    GetOptions($opt, 'type=s', 'iter=n', 'debug');

    my $outputPath = $ARGV[0];
    my $model = LanguageModel::Zerogram->new({ type => $opt->{type}, debug => $opt->{debug} });

    my $istream;
    if ($ARGV[1]) {
	$istream = IO::File->new($ARGV[1]);
    } else {
	$istream = IO::File->new->fdopen(fileno(STDIN), 'r');
    }
    $istream->binmode(':utf8');
    my ($buffer, $id2word) = &readStream($istream, $model);
    foreach my $id (shuffle($buffer)) {
	my $name = $id2word->[$id];
	$model->addWord($name);
    }

    for (my $iter = 1; $iter <= $opt->{iter}; $iter++) {
	printf STDERR ("iter %d\n", $iter) if ($opt->{debug});
	foreach my $id (shuffle($buffer)) {
	    my $name = $id2word->[$id];
	    $model->removeWord($name);
	    $model->addWord($name);
	}
	my $alpha1 = $model->{bigram}->{alpha};
	my $alpha0 = $model->{unigram}->{alpha};
	my $theta1 = $model->{bigram}->{theta};
	my $theta0 = $model->{unigram}->{theta};
	my $d1 = $model->{bigram}->{d};
	my $d0 = $model->{unigram}->{d};
	$model->sampleHyper;

	if ($opt->{debug}) {
	    printf STDERR ("ll: %f\n", $model->calcLogProb);
	    if ($opt->{type} eq 'Dirichlet') {
		printf STDERR ("alpha1 %f -> %f\n", $alpha1, $model->{bigram}->{alpha});
		printf STDERR ("alpha0 %f -> %f\n", $alpha0, $model->{unigram}->{alpha});
	    } else {
		printf STDERR ("theta1 %f -> %f\n", $theta1, $model->{bigram}->{theta});
		printf STDERR ("theta0 %f -> %f\n", $theta0, $model->{unigram}->{theta});
		printf STDERR ("d1 %f -> %f\n", $d1, $model->{bigram}->{d});
		printf STDERR ("d0 %f -> %f\n", $d0, $model->{unigram}->{d});
	    }
	}
    }

    printf STDERR ("save model to %s...\n", $outputPath) if ($opt->{debug});
    use IO::File;
    my $f = IO::File->new($outputPath, 'w') or die;
    $f->binmode(':utf8');
    $model->dump($f);
    $f->close;
}

sub readStream {
    my ($istream, $model) = @_;

    my ($buffer, $id2word, $word2id) = ([], [], {});
    my $counter = 0;
    while ((my $line = $istream->getline)) {
	no warnings qw/uninitialized/;

	chomp($line);
	my ($name, $freq) = split(/\t/, $line, 2);
	my $id;
	if (defined($word2id->{$name})) {
	    $id = $word2id->{$name};
	} else {
	    $id = $word2id->{$name} = $counter;
	    push(@$id2word, $name);
	    $counter++;
	}
	&add2Buffer($buffer, $id, $freq);
    }
    return ($buffer, $id2word);
}

sub add2Buffer {
    my ($buffer, $id, $freq) = @_;

    for (1 .. $freq) {
	push(@$buffer, $id);
    }
}

unless (caller) {
    &main;
}

1;
