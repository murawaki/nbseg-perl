package Text::Sentence::FrequencyBased;

use strict;
use warnings;
use utf8;

use base qw /Text::Sentence/;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	freq => shift,
	list => undef,
    };
    bless($self, $class);
    return $self;
}

sub freq {
    my ($self, $val) = @_;
    if (defined($val)) {
	return $self->{freq} = $val;
    } else {
	return $self->{freq};
    }
}

1;
