package LanguageModel::Table::ZeroAdaptor;
#
# table interface for a zerogram model
#
use strict;
use warnings;
use utf8;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	zerogram => shift,
	isStatic => shift,
    };
    bless($self, $class);
    return $self;
}

sub add {
    my ($self, $word) = @_;
    unless ($self->{isStatic}) {
	$self->{zerogram}->addWord($word);
    }
}
sub remove {
    my ($self, $word) = @_;
    unless ($self->{isStatic}) {
	$self->{zerogram}->removeWord($word);
    }
}
sub draw {
    my ($self) = @_;
    die('draw not implemented yet\n');
}

sub prob {
    my ($self, $name) = @_;
    return $self->{zerogram}->getProb($name);
}

sub logProb {
    my ($self) = @_;
    return $self->{zerogram}->calcLogProb;
}

sub sampleHyper {
    my ($self) = @_;
    $self->{zerogram}->sampleHyper;
}

1;
