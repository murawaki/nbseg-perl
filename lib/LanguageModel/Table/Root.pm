package LanguageModel::Table::Root;
#
# dummy table for root
#
use strict;
use warnings;
use utf8;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	prob => shift,
	drawer => shift,
	isStatic => 1,
    };
    bless($self, $class);
    return $self;
}

sub load {
    my ($this, $self, $drawer) = @_;
    bless($self, $this);
    $self->{drawer} = $drawer;
    $self->{isStatic} = 1;
}

sub TO_JSON {
    my ($self) = @_;
    return {
	prob => $self->{prob},
    };
}

sub prob {
    return ($_[0])->{prob};
}

sub draw {
    my ($self, $context) = @_;
    return &{$self->{drawer}}($context);
}

sub set {
    my ($self, $prob) = @_;
    $self->{prob} = $prob;
}

sub add {}
sub remove {}
sub sampleHyper {}

1;
