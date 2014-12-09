package LanguageModel::Zerogram::UniChar;
#
# zerogram based on character unigram model
#
use strict;
use warnings;
use utf8;

use JSON;

our $POOL_LIMIT = 100000;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	opt => shift,
	ps => 0.5,  # segmentation probability
	a => 1,     # add N smoothing
	total => 0,
	unigram => {},
	pool => {},
	poolLimit => $POOL_LIMIT,
	keyCount => 0,
	time => 0,
    };
    $self->{ps} = $self->{opt}->{ps} if (defined($self->{opt}->{ps}));
    $self->{a} = $self->{opt}->{a}   if (defined($self->{opt}->{a}));
    bless($self, $class);

    return $self;
}

sub addWord {
    my ($self, $name) = @_;

    foreach my $c (split(//, $name)) {
	$self->{unigram}->{$c}++;
	$self->{total}++;
    }
}

sub removeWord {
    my ($self, $name) = @_;

    foreach my $c (split(//, $name)) {
	$self->{unigram}->{$c}--;
	delete($self->{unigram}->{$c}) if ($self->{unigram}->{$c} <= 0);
	$self->{total}--;
    }
}

# sub drawWord {
#     my ($self) = @_;
# }

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

    my $total = $self->{total};
    my $a = $self->{a};
    my $K = scalar(keys(%{$self->{unigram}}));
    my $denom = $total + $K * $a;
    my $M = length($name);

    my $prob = $self->{ps} * ((1 - $self->{ps}) ** ($M - 1));
    foreach my $c (split(//, $name)) {
	$prob *= (($self->{unigram}->{$c} || 0) + $a) / $denom;
    }
    return ($prob > 1e-95)? $prob : 1e-95; # check underflow
}

1;
