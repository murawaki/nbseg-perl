#!/bin/env perl
#
# NOTE: aglign reference and test titles (McNemar's test of statistical significance)
#
use strict;
use warnings;
use utf8;

use Getopt::Long;
use IO::File;

use Text::Sentence::FrequencyBased;
use Text::Sentence::Evaluate qw/evalSentence compareTwoSentences getMcNemar/;

binmode(STDIN,  ':utf8');
binmode(STDOUT, ':utf8');
binmode(STDERR, ':utf8');

my $opt = { };

sub main {
    GetOptions($opt, 'ref=s', 'base=s', 'target=s', 'printDiff');
    my $refList = &loadReference($opt->{ref});
    my $baseList = &loadInput($opt->{base}, $refList);
    my $targetList = &loadInput($opt->{target}, $refList);

    my ($countB, $countC) = (0, 0);
    # while ((my ($raw, $refSentence) = each(%$refList))) {
    foreach my $baseSentence (@$baseList) {
	my $targetSentence = shift(@$targetList);
	my $raw = $baseSentence->asText(1);
	my $raw2 = $targetSentence->asText(1);
	next if ($raw ne $raw2);
	my $refSentence = $refList->{$raw};

	my ($a, $b, $c, $d) = &compareTwoSentences($refSentence, $baseSentence, $targetSentence);
	$countB += $b;
	$countC += $c;

	# NOTE: even if $b == $c, base and target may not the same
	if ($opt->{printDiff}) {
	    if ($b > $c) {
		printf("b > c: %s v. %s (%s)\n", $baseSentence->asText, $targetSentence->asText, $refSentence->asText);
	    } elsif ($b < $c) {
		printf("b < c: %s v. %s (%s)\n", $baseSentence->asText, $targetSentence->asText, $refSentence->asText);
	    }
	}
    }
    print("--------------------------------------------------\n");
    if ($countB + $countC > 0) {
	my $sig = &getMcNemar($countB, $countC);
	printf("McNemar's test: %10f (%s v. %s)\n", $sig, $countB, $countC);
    } else {
	printf("skip McNemar's test because the two outputs are the same\n");
    }
}

sub loadReference {
    my ($path) = @_;

    my $sentenceList = {};
    my $file = IO::File->new($path) or die;
    $file->binmode(':utf8');
    while ((my $line = $file->getline)) {
	chomp($line);
	my ($input, $freq) = split(/\t/, $line);
	$freq = 1 unless ($opt->{freq} && defined($freq));

	$input = '+' . $input . '+' unless ($input =~ /^\+.*\+$/);
	my $sentence = Text::Sentence::FrequencyBased->createFromSegmentedText($input);
	$sentence->freq($freq);
	my $raw = $sentence->asText(1);
	$sentenceList->{$raw} = $sentence;
    }
    $file->close;
    return $sentenceList;
}

sub loadInput {
    my ($path, $wordList) = @_;

    my $sentenceList = [];
    my $file = IO::File->new($path) or die;
    $file->binmode(':utf8');
    while ((my $line = $file->getline)) {
	chomp($line);
	my ($input, $freq) = split(/\t/, $line);
	$freq = 1 unless ($opt->{freq} && defined($freq));
	# early filtering for speed-up
	if (defined($wordList)) {
	    my $raw = $input;
	    $raw =~ s/[\+\-]//g;
	    next unless (defined($wordList->{$raw}));
	}

	$input = '+' . $input . '+' unless ($input =~ /^\+.*\+$/);
	my $sentence = Text::Sentence::FrequencyBased->createFromSegmentedText($input);
	$sentence->freq($freq);
	push(@$sentenceList, $sentence);
    }
    $file->close;
    return $sentenceList;
}

unless (caller) {
    &main;
}

1;
