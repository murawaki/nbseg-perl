package Text::Sentence::Evaluate;

use strict;
use warnings;
use utf8;

use constant {
    B => 0,
    I => 1,
};

our (@ISA, @EXPORT_OK);
BEGIN {
    require Exporter;
    @ISA = qw(Exporter);
    @EXPORT_OK = qw(evalSentence getMcNemar compareTwoSentences getBILabels);
}

use Math::Cephes qw/chdtr/;

sub evalSentence {
    my ($refSentence, $inputSentence) = @_;
    my $refIter = $refSentence->nodeIterator;
    my $inputIter = $inputSentence->nodeIterator;
    my ($correctCount, $refCount, $inputCount) = (0, 0, 0);
    my ($refPos, $inputPos) = (0, 0);
    my ($refNode, $inputNode);
  outer:
    while ((defined ($refNode = $refIter->next))) {
	my $refWord = $refSentence->word($refNode);
	$refCount++;
	while ($refPos > $inputPos) {
	    $inputNode = $inputIter->next;
	    my $inputWord = $inputSentence->word($inputNode);
	    last outer unless (defined($inputNode));
	    $inputCount++;
	    $inputPos += length($inputWord);
	}
	if ($refPos == $inputPos) {
	    $inputNode = $inputIter->next;
	    last unless (defined($inputNode));
	    my $inputWord = $inputSentence->word($inputNode);
	    $inputCount++;
	    $correctCount++ if ($inputWord eq $refWord);
	    $inputPos += length($inputWord);
	}
	$refPos += length($refWord);
    }
    if (defined($refNode)) {
	while ((defined (my $refNode = $refIter->next))) {
	    $refCount++;
	}
    }
    if (defined($inputNode)) {
	while ((defined (my $inputNode = $inputIter->next))) {
	    $inputCount++;
	}
    }
    return ($correctCount, $refCount, $inputCount);
}

sub getMcNemar {
    my ($countB, $countC) = @_;
    return 1.0 - chdtr(1.0, (($countB - $countC) ** 2) / ($countB + $countC));
}

# McNemar's paired test
# sentences are first converted into character-based B/I labels,
# and then labeling disagreements are examined
# see [Kudo+ 2004]
sub compareTwoSentences {
    my ($refSentence, $inputSentence1, $inputSentence2) = @_;
    my $refLabels = &getBILabels($refSentence);
    my $inputLabels1 = &getBILabels($inputSentence1);
    my $inputLabels2 = &getBILabels($inputSentence2);

    # $a: both correct 
    # $b: 1 correct, 2 wrong
    # $c: 1 wrong, 2 correct
    # $d: both wrong
    my ($a, $b, $c, $d) = (0, 0, 0, 0);
    foreach my $rl (@$refLabels) {
	my $il1 = shift(@$inputLabels1);
	my $il2 = shift(@$inputLabels2);
	if ($rl == $il1) {
	    if ($il1 == $il2) {
		$a++;
	    } else {
		$b++;
	    }
	} else {
	    if ($rl == $il2) {
		$c++;
	    } else {
		$d++;
	    }
	}
    }
    return ($a, $b, $c, $d);
}

sub getBILabels {
    my ($sentence) = @_;
    my $rv = [];
    my $iter = $sentence->nodeIterator;
    while ((my $node = $iter->next)) {
	my $l = length($sentence->word($node));
	push(@$rv, B);
	for (1 .. $l - 1) {
	    push(@$rv, I);
	}
    }
    return $rv;
}

1;
