package LanguageModel::Table::PitmanYorProcess;
#
# a group of Pitman-Yor processes
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
our ($THETA_A, $THETA_B) = (0.01, 0.01); # G(a,b), the prior for theta
our ($D_A, $D_B) = (0.1, 0.1);           # B(a,b), the piior for d

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {
	backTable => shift,
	contextReducer => shift,  # ref to a subroutine
	isStatic => 0,
	restaurantList => {},
	theta => 5, d => 0.8, # default values for hyperparameters
	thetaA => $THETA_A, thetaB => $THETA_B,
	dA => $D_A, dB => $D_B,
    };
    bless($self, $class);
    return $self;
}

sub load {
    my ($this, $self, $backTable, $reducer) = @_;

    $self->{backTable} = $backTable;
    $self->{contextReducer} = $reducer;
    ($self->{dA}, $self->{dB}) = ($D_A, $D_B);
    ($self->{thetaA}, $self->{thetaB}) = ($THETA_A, $THETA_B);
    bless($self, $this);
}

sub TO_JSON {
    my ($self) = @_;

    # drop 'backTable' to avoid recursive link
    return {
	isStatic => 0,
	restaurantList => $self->{restaurantList},
	theta => $self->{theta},
	d => $self->{d},
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
    my ($d, $theta) = ($self->{d}, $self->{theta});
    my $tableCount = 0;
    my $histogram = $restaurant->{histogram}->{$current};
    if (defined($histogram)) {
	$tableCount = scalar(@$histogram);
    }

    my $denom = $theta + $restaurant->{marginal};
    my $num1 = ($restaurant->{count}->{$current} || 0) - $d * $tableCount;
    my $num2 = ($theta + $d * $restaurant->{tableMarginal})
	* $self->{backTable}->prob(&{$self->{contextReducer}}($current, $context));
    my $prob = ($num1 + $num2) / $denom;
    return ($prob > 1e-95)? $prob : 1e-95; # check underflow
}

sub draw {
    my ($self, $context) = @_;

    my $d = $self->{d};
    my $restaurant = $self->getRestaurant($context);
    my $pShare = $restaurant->{marginal} - $d * $restaurant->{tableMarginal};
    my $pNew = $self->{theta} + $d * $restaurant->{tableMarginal};
    my $r = gsl_ran_flat($rng->raw, 0, $pShare + $pNew);
    if ($r < $pShare) {
	my $r = gsl_ran_flat($rng->raw, 0, $pShare);
	my $current;
	my $countList = $restaurant->{count};
	my @currentList = keys(%$countList);
	foreach $current (@currentList) {
	    $r -= $countList->{$current} - $d;
	    return $current if ($r <= 0);
	}
	return $currentList[-1];
    } else {
	my (undef, $context2) = &{$self->{contextReducer}}(undef, $context);
	return $self->{backTable}->draw($context2);
    }
}

sub add {
    my ($self, $current, $context) = @_;
    my ($current2, $context2) = &{$self->{contextReducer}}($current, $context);

    my $restaurant = $self->getRestaurant($context, 1); # create when undefined
    my $newTable = 1;
    my $histogram = $restaurant->{histogram}->{$current};
    if (defined($histogram)) {
	my $d = $self->{d};
	my $pShare = $restaurant->{count}->{$current} - $d * scalar(@$histogram);
	my $pNew = ($self->{theta} + $d * ($restaurant->{tableMarginal} || 0))
	    * $self->{backTable}->prob($current2, $context2);
	my $r = gsl_ran_flat($rng->raw, 0, $pShare + $pNew);
	if ($r < $pShare) {
	    $newTable = 0;
	    my $added = 0;
	    for (my $i = 0, my $l = scalar(@$histogram); $i < $l; $i++) {
		my $c = $histogram->[$i];
		$r -= $c - $d;
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
	# for $pNew: go through
    }
    if ($newTable) {
	push(@{$restaurant->{histogram}->{$current}}, 1);
	$restaurant->{tableMarginal}++;
	$self->{backTable}->add($current2, $context2);
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
    my $theta = $self->{theta};
    my $d = $self->{d};
    my $logD = log($d);
    my $lgamTheta = lgam($theta);
    my $lgamThetaD = lgam($theta / $d);
    my $lgamOneMinusD = lgam(1 - $d);

    while ((my ($context, $restaurant) = each(%{$self->{restaurantList}}))) {
	# collecting denominators
	$ll += $lgamTheta - lgam($theta + $restaurant->{marginal});

	# prob. of creating new tables
	my $tableCount = $restaurant->{tableMarginal};
	# generalized factorial
	$ll += $tableCount * $logD + lgam($theta / $d + $tableCount) - $lgamThetaD;

	# prob. of sitting at existing tables
	my $histogramList = $restaurant->{histogram};
	while ((my ($current, $histogram) = each(%$histogramList))) {
	    foreach my $count (@$histogram) {
		next if ($count <= 1);
		$ll += lgam($count - $d) - $lgamOneMinusD;
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

sub theta {
    my ($self, $newVal) = @_;
    if (defined($newVal)) {
	return $self->{theta} = $newVal;
    } else {
	return $self->{theta};
    }
}

sub d {
    my ($self, $newVal) = @_;
    if (defined($newVal)) {
	return $self->{d} = $newVal;
    } else {
	return $self->{d};
    }
}

sub sampleHyper {
    my ($self) = @_;

    my ($dOld, $thetaOld) = ($self->{d}, $self->{theta});
    my ($dList, $thetaList) = ([], []);
    for (1 .. $HP_DRAW_COUNT) {
	($dOld, $thetaOld) = $self->sampleHyperSub($dOld, $thetaOld);
	push(@$dList, $dOld);
	push(@$thetaList, $thetaOld);
    }
    $self->{d} = sum(@$dList) / $HP_DRAW_COUNT;
    $self->{theta} = sum(@$thetaList) / $HP_DRAW_COUNT;
}


sub sampleHyperSub {
    my ($self, $dOld, $thetaOld) = @_;

    # hyperparameter sampling based on:
    #   Teh+ 2006: A Bayesian Interpretation of Interpolated Kneser-Ney

    my ($x, $y, $Y, $z) = (0, 0, 0, 0);
    while ((my ($context, $restaurant) = each(%{$self->{restaurantList}}))) {
	my $marginalCount = $restaurant->{marginal};
	if ($marginalCount >= 2) {
	    $x += log(gsl_ran_beta($rng->raw, $thetaOld + 1, $marginalCount - 1));
	}

	my $tableMarginal = $restaurant->{tableMarginal};
    	for (my $i = 1; $i < $tableMarginal; $i++) {
    	    $y += (gsl_ran_flat($rng->raw, 0, $thetaOld + $dOld * $i) < $thetaOld)? 1 : 0;
	    $Y++;
    	}

	my $histogramList = $restaurant->{histogram};
	while ((my ($current, $histogram) = each(%$histogramList))) {
	    foreach my $c (@$histogram) {
		for (my $j = 1; $j < $c; $j++) {
		    $z += (gsl_ran_flat($rng->raw, 0, $j - $dOld) < $j - 1)? 0 : 1; # reverse
		}
	    }
	}
    }
    my $dNew = gsl_ran_beta($rng->raw, $self->{dA} + $Y - $y, $self->{dB} + $z);
    my $thetaNew = gsl_ran_gamma($rng->raw, $self->{thetaA} + $y, 1 / ($self->{thetaB} - $x));
    return ($dNew, $thetaNew);
}

1;
