#!/bin/env perl
#
# load a zerogram (character n-gram) model
# draws a word/returns the probability to a given word
# 
# usage: program MODEL
#
use strict;
use warnings;
use utf8;

use IO::File;
use LanguageModel::Zerogram;
use Dumpvalue;

binmode(STDIN,  ':utf8');
binmode(STDOUT, ':utf8');
binmode(STDERR, ':utf8');

my $modelPath = $ARGV[0];
my $model;
{
    my $f = IO::File->new($modelPath);
    $f->binmode(':utf8');
    $model = LanguageModel::Zerogram->load($f);
    $f->close;
}

my $iostream;
if ($ARGV[1]) {
    $iostream = IO::File->new($ARGV[1]);
} else {
    $iostream = IO::File->new->fdopen(fileno(STDIN), 'r');
}
$iostream->binmode(':utf8');
while ((my $line = $iostream->getline)) {
    no warnings qw/uninitialized/;

    chomp($line);
    if ($line eq 'D') {
	printf("%s\n", $model->drawWord);
    } else {
	printf("%s\t%.30f\n", $line, $model->getProb($line));
    }
}

1;
