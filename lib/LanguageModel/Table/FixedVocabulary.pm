package LanguageModel::Table::FixedVocabulary;
#
# Fixed vocabulary modeled with a Dirichlet distribution
#
use strict;
use warnings;
use utf8;

use Carp::Assert;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	vocSize => undef,
	count => {},         # word -> count
	marginal => 0,
	alpha => 1.0,
    };
    bless($self, $class);
    return $self;
}

sub vocSize {
    my ($self) = @_;

    if (defined($self->{vocSize})) {
	return $self->{vocSize};
    } else {
	return $self->{vocSize} = keys(%{$self->{count}});
    }
}

sub prob {
    my ($self, $name) = @_;

    return (($self->{count}->{$name} || 0) + $self->{alpha})
	/ ($self->{marginal} + $self->{alpha} * $self->vocSize);
}

sub add {
    my ($self, $name) = @_;

    {
	no warnings qw/uninitialized/;
	$self->{count}->{$name}++;
    }
    $self->{marginal}++;
    undef($self->{vocSize});
}

sub remove {
    my ($self, $name) = @_;

    $self->{count}->{$name}--;
    assert($self->{count}->{$name} >= 0, sprintf("negative count: %d", $name));
    if ($self->{count}->{$name} == 0) {
	delete($self->{count}->{$name});
    }
    $self->{marginal}--;
    assert($self->{marginal} >= 0, "negative marginal count");
    undef($self->{vocSize});
}

1;
