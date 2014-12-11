package Text::Sentence::MergerSplitter;
#
# hold a pointer to a node
# merge/split operations
#
# Note: In the current implementation, the merge/split operation in
# the same sentence does not affect the current node of the
# char. iterator.
#
use strict;
use warnings;
use utf8;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	sentence => shift,
	node => shift,     # pointer to the current node
	wpos => shift,     # character-wise pos. of the current node in the sentence
    };
    bless($self, $class);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    undef($self->{sentence});
    undef($self->{node});
}

sub wpos {
    my ($self) = @_;
    return $self->{wpos};
}

# merge the current node with the right node
# returns the new node
sub merge {
    my ($self, $newWord) = @_;

    my $leftNode = $self->{node};
    my $rightNode = $self->{node}->[Text::Sentence->NEXT];
    $self->{node} = $self->{sentence}->add($rightNode, $self->{sentence}->segType($rightNode), $newWord);
    $self->{sentence}->delete($leftNode);
    $self->{sentence}->delete($rightNode);
    return $self->{node};
}

sub split {
    my ($self, $newLeft, $newRight) = @_;

    my $currentNode = $self->{node};
    $self->{node} = $self->{sentence}->add($currentNode, Text::Sentence->BOUNDARY, $newLeft);
    $self->{sentence}->add($self->{node}, $self->{sentence}->segType($currentNode), $newRight);
    $self->{sentence}->delete($currentNode);
    return $self->{node};
}

1;


package Text::Sentence::CharacterIterator;
#
# NOTE: the current node is the left one when examining a boundary
#
use strict;
use warnings;
use utf8;
use base qw/Text::Sentence::MergerSplitter/;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	sentence => shift,
	nodePos => 0,
	wpos => 0,
    };
    $self->{node} = $self->{sentence}->{list};
    bless($self, $class);
    return $self;
}

sub next {
    my ($self) = @_;

    my $sentence = $self->{sentence};
    my $word = $sentence->word($self->{node});
    my $len = length($word);
    if (++$self->{nodePos} < $len) {
	# inside a node
	return {
	    node => $self->{node},
	    segType => Text::Sentence->NONE_BOUNDARY,
	    left => substr($word, 0, $self->{nodePos}),
	    right => substr($word, $self->{nodePos}),
	    unsegmented => $word,
	    isRightmost => (defined($self->{node}->[Text::Sentence->NEXT]))? 0 : 1,
	};
    } elsif ($self->{nodePos} > $len) {
	# the next node
	$self->{nodePos} = 0;
	$self->{wpos} += $len;
	$self->{node} = $self->{node}->[Text::Sentence->NEXT];
	unless (defined($self->{node})) {
	    print STDERR ("sentence iteration out of range\n");
	    return undef;
	}
	return $self->next;
    } else {
	# boundary
	my $rightNode = $self->{node}->[Text::Sentence->NEXT];
	return undef unless (defined($rightNode));
	if ($sentence->segType($self->{node}) == Text::Sentence->FIXED_BOUNDARY) {
	    return $self->next;
	}
	my $right = $sentence->word($rightNode);
	return {
	    node => $self->{node},
	    segType => Text::Sentence->BOUNDARY,
	    left => $word,
	    right => $right,
	    unsegmented => $word . $right,
	    isRightmost => (defined($rightNode->[Text::Sentence->NEXT]))? 0 : 1,
	};
    }
}

# merge the current node with the right node
# returns the new node
sub merge {
    my ($self, $newWord) = @_;
    my $oldLeft = $self->{sentence}->word($self->{node});
    $self->{nodePos} = length($oldLeft);
    return $self->SUPER::merge($newWord);
}
# no change needed for split

1;


package Text::Sentence::NodeIterator;

use strict;
use warnings;
use utf8;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	sentence => shift,
	wpos => 0,
	current => undef,
    };
    bless($self, $class);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    undef($self->{sentence});
}

sub next {
    my ($self) = @_;
    if (defined($self->{current})) {
	$self->{wpos} += length($self->{sentence}->word($self->{current}));
	return $self->{current} = $self->{current}->[Text::Sentence->NEXT];
    } else {
	return ($self->{current} = $self->{sentence}->{list});
    }
}

sub wpos {
    my ($self) = @_;
    return $self->{wpos};
}

1;

package Text::Sentence;

use strict;
use warnings;
use utf8;

use base qw /Class::Data::Inheritable/;

use Math::GSL::RNG qw/gsl_rng_uniform/;
use LanguageModel::Util qw/$rng/;

__PACKAGE__->mk_classdata(NONE_BOUNDARY => 0);
__PACKAGE__->mk_classdata(BOUNDARY => 1);
__PACKAGE__->mk_classdata(FIXED_BOUNDARY => 2);

__PACKAGE__->mk_classdata(PREV => 0);
__PACKAGE__->mk_classdata(NEXT => 1);
__PACKAGE__->mk_classdata(SEG_TYPE => 2); # of the right boundary
__PACKAGE__->mk_classdata(WORD => 3);

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	list => undef,
    };
    bless($self, $class);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    $self->deconstructList;
    delete($self->{list});
}

sub sid {
    my ($self, $val) = @_;
    if (defined($val)) {
	return $self->{sid} = $val;
    } else {
	return $self->{sid};
    }
}

sub segType {
    my ($self, $node, $val) = @_;
    if (defined($val)) {
	return $node->[$self->SEG_TYPE] = $val;
    } else {
	return $node->[$self->SEG_TYPE];
    }
}

sub word {
    my ($self, $node, $val) = @_;
    if (defined($val)) {
	return $node->[$self->WORD] = $val;
    } else {
	return $node->[$self->WORD];
    }
}

# left_node, seg_type, string
sub add {
    my ($self, $a, $seg, $word) = @_;
    my $node = [undef, undef, $seg, $word];
    if (defined($a)) {
	my $b = $a->[$self->NEXT];
	$a->[$self->NEXT] = $node;
	$node->[$self->PREV] = $a;
	$node->[$self->NEXT] = $b;
	$b->[$self->PREV] = $node if ($b);
    } else {
	$self->{list} = $node;
    }
    return $node;
}

# node_to_be_deleted
sub delete {
    my ($self, $a) = @_;
    my $p = $a->[$self->PREV];
    my $n = $a->[$self->NEXT];
    $a->[$self->PREV] = $a->[$self->NEXT] = undef;
    if (defined($p)) {
	$p->[$self->NEXT] = $n;	
	$n->[$self->PREV] = $p if ($n);
    } else {
	if (defined($n)) {
	    $self->{list} = $n;
	    $n->[$self->PREV] = undef;
	} else {
	    $self->{list} = undef;
	}
    }
}

sub replaceFully {
    my ($self, $wordList) = @_;

    $self->deconstructList;    
    my $node = undef;
    foreach my $tmp (@$wordList) {
	my ($word, $segType) = @$tmp;
	$node = $self->add($node, $segType, $word);
    }
}

sub createFromSegmentedText {
    my $self = (shift)->new;
    my $segText = shift;
    my $node = undef;
    foreach my $seg (split(/(?<!\\)\+/, $segText)) {
	$seg =~ s/\\\+/\+/g;
	foreach my $seg2 (split(/(?<!\\)\-/, $seg)) {
	    $seg2 =~ s/\?//g;
	    $seg2 =~ s/\\\-/\-/g;
	    $node = $self->add($node, ref($self)->BOUNDARY, $seg2);
	}
	$self->segType($node, ref($self)->FIXED_BOUNDARY);
    }
    return $self;
}

sub nodeIterator {
    my ($self) = @_;
    return Text::Sentence::NodeIterator->new($self);
}

sub charIterator {
    my ($self) = @_;
    return Text::Sentence::CharacterIterator->new($self);
}

sub getNodeHolder {
    my ($self, $node, $wpos) = @_;
    return Text::Sentence::MergerSplitter->new($self, $node, $wpos);
}

sub length {
    my ($self) = @_;

    if (defined(my $l = $self->{length})) {
	return $l;
    }
    return ($self->{length} = length($self->asText(1)));
}

sub deconstructList {
    my ($self) = @_;
    my $node = $self->{list};
    while (1) {
	my $next = $node->[$self->NEXT];
	undef($node->[$self->PREV]);
	undef($node->[$self->NEXT]);
	last unless (defined($next));
	$node = $next;
    }
}

# randomize segmentation
# $eta: boundary probability
#
# NOTE: use only at the time of initialization
#       otherwise iterators get broken
sub randomize {
    my ($self, $eta) = @_;

    my $sentence = ref($self)->new;
    my $cnode = undef;
    my $iter = $self->nodeIterator;
    my $remnant = '';
    while ((my $node = $iter->next)) {
	my $word = $node->[$self->WORD];
	my $segType = $node->[$self->SEG_TYPE];
	my @charList = split(//, $word);
	my $last;
	if ($segType == $self->FIXED_BOUNDARY) {
	    $last = pop(@charList);
	}
	foreach my $c (@charList) {
	    if (gsl_rng_uniform($rng->raw) <= $eta) {
		$cnode = $sentence->add($cnode, $self->BOUNDARY, $remnant . $c);
		$remnant = '';
	    } else {
		$remnant .= $c;
	    }
	}
	if (defined($last)) {
	    $cnode = $sentence->add($cnode, $self->FIXED_BOUNDARY, $remnant . $last);
	    $remnant = '';
	}
    }
    if (CORE::length($remnant) > 0) {
	$cnode = $sentence->add($cnode, $self->FIXED_BOUNDARY, $remnant);
    }
    # HACK
    $self->deconstructList;
    $self->{list} = $sentence->{list};
    undef($sentence->{list});
}

sub asText {
    my ($self, $plain) = @_;

    my $FIXBND = '+';
    my $BND = '-';
    if ($plain) {
	$FIXBND = $BND = '';
    }
    my $rv = $FIXBND;
    my $iter = $self->nodeIterator;
    while ((my $node = $iter->next)) {
	my $word = $node->[$self->WORD];
	unless ($plain) {
	    $word =~ s/\-/\\-/;
	    $word =~ s/\+/\\+/;
	}
	$rv .= $word;
	$rv .= ($node->[$self->SEG_TYPE] == ref($self)->BOUNDARY)? $BND : $FIXBND;
    }
    return $rv;
}

# for debug
sub checkConsistency {
    my ($self) = @_;

    my $node = $self->{sentence}->{list};
    if (defined($node->[Text::Sentence->PREV])) {
	printf STDERR ("first node has a previous node: %s", $node->[Text::Sentence->WORD]);
    }
    my $prev = $node;
    while ((my $node = $node->[Text::Sentence->NEXT])) {
	unless ($prev->[Text::Sentence->NEXT] == $node) {
	    printf STDERR ("inconsistent %s -> %s", $prev->[Text::Sentence->WORD], $node->[Text::Sentence->WORD]);
	}
	unless ($node->[Text::Sentence->PREV] == $prev) {
	    printf STDERR ("inconsistent %s <- %s", $prev->[Text::Sentence->WORD], $node->[Text::Sentence->WORD]);
	}
	$prev = $node;
    }
    # TODO: check cycle?
}

1;
