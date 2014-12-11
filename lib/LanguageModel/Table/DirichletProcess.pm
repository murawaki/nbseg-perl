package LanguageModel::Table::DirichletProcess;
#
# a group of Dirichlet processes
#
# NOTE: explicit table tracking is unnecessary for the unigram model
#
use strict;
use warnings;
use utf8;

use Carp::Assert;
use List::Util qw/sum/;
use Math::Cephes qw/lgam psi/;
use Math::GSL::Randist qw/gsl_ran_flat gsl_ran_beta gsl_ran_gamma/;

use LanguageModel::Util qw/$rng/;

our $HP_DRAW_COUNT = 5;
our ($ALPHA_A, $ALPHA_B) = (0.01, 0.01); # G(a,b), the prior for the hyperparameter

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	backTable => shift,
	contextReducer => shift,  # ref to a subroutine
	isStatic => 0,
	restaurantList => {},
	alpha => 1,
	a => $ALPHA_A, b => $ALPHA_B,
    };
    bless($self, $class);
    return $self;
}

sub load {
    my ($this, $self, $backTable, $reducer) = @_;

    $self->{backTable} = $backTable;
    $self->{contextReducer} = $reducer;
    ($self->{a}, $self->{b}) = ($ALPHA_A, $ALPHA_B);
    bless($self, $this);
}

sub TO_JSON {
    my ($self) = @_;

    # drop 'backTable' to avoid recursive link
    return {
	isStatic => 0,
	restaurantList => $self->{restaurantList},
	alpha => $self->{alpha},
    };
}

sub getRestaurant {
    my ($self, $context, $doCreate) = @_;
    $doCreate = 0 unless (defined($doCreate));

    my $restaurant = $self->{restaurantList}->{$context};
    unless (defined($restaurant)) {
	$restaurant = {
	    histogram => {},
	    count => {},
	    marginal => 0,
	    tableMarginal => 0,
	};
	$self->{restaurantList}->{$context} = $restaurant if ($doCreate);
    }
    return $restaurant;
}

sub prob {
    my ($self, $current, $context) = @_;

    my $restaurant = $self->getRestaurant($context);
    my $alpha = $self->{alpha};
    my $prob = (($restaurant->{count}->{$current} || 0) + $alpha * $self->{backTable}->prob(&{$self->{contextReducer}}($current, $context)))
	/ ($restaurant->{marginal} + $alpha);
    return ($prob > 1e-95)? $prob : 1e-95; # check underflow
}

sub draw {
    my ($self, $context) = @_;

    my $restaurant = $self->getRestaurant($context);
    my $pShare = $restaurant->{marginal};
    my $pNew = $self->{alpha};
    if ($pShare <= 0 || gsl_ran_flat($rng->raw, 0, $pShare + $pNew) < $pNew) {
	my (undef, $context2) = &{$self->{contextReducer}}(undef, $context);
	return $self->{backTable}->draw($context2);
    } else {
	my $r = gsl_ran_flat($rng->raw, 0, $pShare);
	my $current;
	my $countList = $restaurant->{count};
	my @currentList = keys(%$countList);
	foreach $current (@currentList) {
	    $r -= $countList->{$current};
	    return $current if ($r <= 0);
	}
	return $currentList[-1];
    }
}

sub add {
    my ($self, $current, $context) = @_;
    my ($current2, $context2) = &{$self->{contextReducer}}($current, $context);

    my $restaurant = $self->getRestaurant($context, 1); # create when undefined
    # update histogram
    my $pShare = ($restaurant->{count}->{$current} || 0);
    my $pNew = $self->{alpha} * $self->{backTable}->prob($current2, $context2);
    if ($pShare <= 0 || gsl_ran_flat($rng->raw, 0, $pShare + $pNew) < $pNew) {
	push(@{$restaurant->{histogram}->{$current}}, 1);
	$restaurant->{tableMarginal}++;
	$self->{backTable}->add($current2, $context2);
    } else {
	my $r = gsl_ran_flat($rng->raw, 0, $pShare);
	my $added = 0;
	my $histogram = $restaurant->{histogram}->{$current};
	for (my $i = 0, my $l = scalar(@$histogram); $i < $l; $i++) {
	    my $c = $histogram->[$i];
	    $r -= $c;
	    if ($r <= 0) {
		$histogram->[$i]++;
		$added = 1;
		last;
	    }
	}
	unless ($added) { # fail safe
	    $histogram->[-1]++;
	}
    }
    # update count
    $restaurant->{count}->{$current}++;
    # update marginal
    $restaurant->{marginal}++;
}

sub remove {
    my ($self, $current, $context) = @_;

    my $restaurant = $self->getRestaurant($context);
    my $count = $restaurant->{count}->{$current};
    assert(defined($count) && $count > 0, sprintf("[removal] non-existent word: %s:%s", $context, $current)) if (DEBUG);

    my $r = gsl_ran_flat($rng->raw, 0, $count);
    my $deleted = 0; my $emptyIndex = -1;
    my $histogram = $restaurant->{histogram}->{$current};
    for (my $i = 0, my $l = scalar(@$histogram); $i < $l; $i++) {
	my $c = $histogram->[$i];
	$r -= $c;
	if ($r <= 0) {
	    if (--$histogram->[$i] <= 0) {
		$emptyIndex = $i;
	    }
	    $deleted = 1;
	    last;
	}
    }
    unless ($deleted) { # fail safe
	$histogram->[-1]--;
	if ($histogram->[-1] <= 0) {
	    $emptyIndex = scalar(@$histogram) - 1;
	}
    }
    my $isEmpty = 0;
    if ($emptyIndex >= 0) {
	splice(@$histogram, $emptyIndex, 1);
	my $tableMarginal = --$restaurant->{tableMarginal};
	if ($tableMarginal <= 0) {
	    $isEmpty = 1;
	}
	$self->{backTable}->remove(&{$self->{contextReducer}}($current, $context));
    }
    if (scalar(@$histogram) <= 0) {
	delete($restaurant->{histogram}->{$current});
	assert(scalar(keys(%{$restaurant->{histogram}})) > 0 || $isEmpty) if (DEBUG);
    }

    # update count
    $count = --$restaurant->{count}->{$current};
    if ($count <= 0) {
	delete($restaurant->{count}->{$current});
	assert(scalar(keys(%{$restaurant->{count}})) > 0 || $isEmpty) if (DEBUG);
    }

   # update marginal
    my $marginal = --$restaurant->{marginal};
    assert($marginal > 0 || $isEmpty) if (DEBUG);
    if ($isEmpty) {
	delete($self->{restaurantList}->{$context});
    }
}

sub logProb {
    my ($self) = @_;
    my $ll = 0;
    my $alpha = $self->{alpha};
    my $logAlpha = log($alpha);
    my $lgamAlpha = lgam($alpha);

    while ((my ($context, $restaurant) = each(%{$self->{restaurantList}}))) {
	# collecting denominators
	$ll += $lgamAlpha - lgam($alpha + $restaurant->{marginal});

	# prob. of creating new tables
	$ll += $restaurant->{tableMarginal} * $logAlpha;

	# prob. of sitting at existing tables
	my $histogramList = $restaurant->{histogram};
	while ((my ($current, $histogram) = each(%$histogramList))) {
	    foreach my $count (@$histogram) {
		next if ($count <= 1);
		$ll += lgam($count);  # (count - 1)!
	    }
	}
    }
    # prob. of backoff generation
    if ($self->{backTable}->{isStatic}) {
	while ((my ($context, $restaurant) = each(%{$self->{restaurantList}}))) {
	    my $histogramList = $restaurant->{histogram};
	    while ((my ($current, $histogram) = each(%$histogramList))) {
		$ll += scalar(@$histogram) * log($self->{backTable}->prob(&{$self->{contextReducer}}($current, $context)));
	    }
	}
	return $ll;
    } else {
	return $ll + $self->{backTable}->logProb;
    }
}

sub alpha {
    my ($self, $newVal) = @_;
    if (defined($newVal)) {
	return $self->{alpha} = $newVal;
    } else {
	return $self->{alpha};
    }
}

sub sampleHyper {
    my ($self) = @_;

    # MCMC estimation
    my $alphaOld = $self->{alpha};
    my $alphaList = [];
    for (1 .. $HP_DRAW_COUNT) {
	$alphaOld = $self->drawAlpha($alphaOld);
	push(@$alphaList, $alphaOld);
    }
    $self->{alpha} = sum(@$alphaList) / $HP_DRAW_COUNT;
}

sub drawAlpha {
    my ($self, $alphaOld, $marginalIndex) = @_;

    # another hyperparameter sampling based on
    #   Teh et al.: Hierarchical Dirichlet Process
    my $a1 = $self->{a};
    my $b1 = $self->{b};

    while ((my ($context, $restaurant) = each(%{$self->{restaurantList}}))) {
	my $T = $restaurant->{tableMarginal};
	my $n = $restaurant->{marginal};
	$a1 += $T - ((gsl_ran_flat($rng->raw, 0, $alphaOld) < $n)? 1 : 0);
	$b1 -= log(gsl_ran_beta($rng->raw, $alphaOld + 1, $n));
    }
    return gsl_ran_gamma($rng->raw, $a1, 1 / $b1);
}

sub drawAlphaFlat {
    my ($self, $alphaOld, $k, $n) = @_;
    # k: total number of tables
    # n: total number of instances

    # hyperparameter sampling for flat (non-hierarchical) model
    #   Michael D. Escobar and Mike West:
    #   Bayesian Density Estimation and Inference Using Mixtures
    my $a1 = $self->{a} + $k;
    my $a2 = $self->{a} + $k - 1;

    my $logEta = log(gsl_ran_beta($rng->raw, $alphaOld + 1, $n));
    my $pi = ($self->{a} + $k - 1) / ($self->{a} + $k + $n * $self->{b} - (1 + $n * $logEta));
    my $b1 = $self->{b} - $logEta;
    return $pi * gsl_ran_gamma($rng->raw, $a1, 1 / $b1)
	+ (1.0 - $pi) * gsl_ran_gamma($rng->raw, $a2, 1 / $b1);
}

1;
