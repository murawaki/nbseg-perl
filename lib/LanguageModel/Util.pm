package LanguageModel::Util;

use strict;
use warnings;
use utf8;

our (@ISA, @EXPORT_OK);
BEGIN {
    require Exporter;
    @ISA = qw(Exporter);
    @EXPORT_OK = qw($rng setSeed shuffle randPartition randPartitionLog logsumexp logsumexpList getFlipProposal sliceSampler1d);
}

use Carp::Assert;
use List::Util qw/min max sum/;
use Math::GSL::RNG qw/$gsl_rng_mt19937 gsl_rng_uniform gsl_rng_uniform_int/;
use Math::GSL::Randist qw/gsl_ran_flat gsl_ran_beta_pdf gsl_ran_ugaussian_pdf/;
use Math::Cephes qw/lgam/;

our $rng = Math::GSL::RNG->new($gsl_rng_mt19937, time ^ $$);

sub setSeed {
    $rng = Math::GSL::RNG->new($gsl_rng_mt19937, $_[0]);
}

# destructive shuffle
sub shuffle {
    return @_ if !@_ || ref $_ [0] eq 'ARRAY' && !@{$_ [0]};
    my $array = @_ == 1 && ref $_ [0] eq 'ARRAY' ? shift : [@_];
    for (my $i = @$array; -- $i;) {
        my $r = gsl_rng_uniform_int($rng->raw, $i);
       ($array->[$i], $array->[$r]) = ($array->[$r], $array->[$i]);
    }
    return wantarray? @$array : $array;
}

sub randPartition {
    my ($probArray) = @_;

    my $sum = 0.0;
    foreach my $p (@$probArray) {
	$sum += $p;
    }
    my $r = gsl_ran_flat($rng->raw, 0, $sum);
    for (my $i = 0, my $l = scalar(@$probArray); $i < $l; $i++) {
	$r -= $probArray->[$i];
	if ($r <= 0) {
	    return $i;
	}
    }
    return scalar(@$probArray) - 1;
}

sub randPartitionLog {
    my ($logProbArray, $anneal) = @_;
    $anneal = 1 unless (defined($anneal));

    my $massList = [];
    my $sum = 0.0;

    my $maxI = 0;
    my $maxVal = $logProbArray->[0];
    for (my $i = 1, my $l = scalar(@$logProbArray); $i < $l; $i++) {
	if ($logProbArray->[$i] > $maxVal) {
	    $maxI = $i;
	    $maxVal = $logProbArray->[$i];
	}
    }

    my $base = int(-1 * $logProbArray->[$maxI]);
    foreach my $l (@$logProbArray) {
	my $n = exp($l + $base);
	$sum += $n;
	push(@$massList, $n);
    }
    if ($anneal != 1) {
	my $sum2 = 0;
	foreach my $mass (@$massList) {
	    my $v = ($mass / $sum) ** $anneal;
	    $mass = $v;
	    $sum2 += $v;
	}
	$sum = $sum2;
    }

    my $rand = gsl_ran_flat($rng->raw, 0, $sum);
    for (my $i = 0, my $l = scalar(@$massList); $i < $l; $i++) {
	$rand -= $massList->[$i];
	if ($rand <= 0) {
	    return $i;
	}
    }
    return scalar(@$massList) - 1;
}

sub logsumexp {
    my ($a, $b) = @_;
    ($a, $b) = ($a < $b)? ($a, $b) : ($b, $a);
    return $b + log(exp($a-$b) + 1.0);
}

sub logsumexpList {
    my ($list) = @_;
    my $c = min(@$list);
    return $c + log(sum(map { exp($_ - $c) } (@$list)));
}

our $PROB_PERMUTATION = 0.1;
our $PROPOSAL_LAMBDA = 0.5;

# proposal distribution
#   discretized mixture of normal and beta distributions
#   the former prefers local poinsts while the latter try to move to the edges
sub getFlipProposal {
    my ($current, $M) = @_;

    # just permutation
    if ($current != 0 && $current != $M) {
	if (gsl_ran_flat($rng->raw, 0, 1) < $PROB_PERMUTATION) {
	    return ($current, 1, 1);
	}
    }

    # prob. mass table of Q(h'|h)
    my $mean = $current / $M;
    my $sigma = sqrt($M / 4);
    my $sum = 0.0;
    my $mass = [];
    my $bpList = []; # beta distribution is symmetric
    foreach my $i (0 .. $M) {
	my $point = $i / $M;
	if ($i == 0) {
	    $point += 0.05 / $M;
	} elsif ($i == $M) {
	    $point -= 0.05 / $M;
	}
	my $bp = gsl_ran_beta_pdf($point, 0.1, 0.1);
	my $np = gsl_ran_ugaussian_pdf(($i - $current) / $sigma);
	my $val = $PROPOSAL_LAMBDA * $bp + (1 - $PROPOSAL_LAMBDA) * $np;
	{
	    no warnings qw/uninitialized/;
	    $bpList->[$i] = $bp;
	    $mass->[$i] = $val;
	}
	if ($i != $current) {
	    $sum += $val;
	}
    }
    my $proposal;
    my $r = gsl_ran_flat($rng->raw, 0, $sum);
    foreach my $i (0 .. $M) {
	next if ($i == $current);
	$r -= $mass->[$i];
	if ($r <= 0) {
	    $proposal = $i;
	    last;
	}
    }
    unless (defined($proposal)) { # fail-safe
	$proposal = ($current == $M)? $M : $M - 1;
    }

    # now calc Q(h|h')
    my $logQF = log($mass->[$proposal] / $sum);
    $sum = 0.0;
    $mean = $proposal / $M;
    foreach my $i (0 .. $M) {
	my $point = $i / $M;
	my $bp = $bpList->[$i];
	my $np = gsl_ran_ugaussian_pdf(($i - $current) / $sigma);
	my $val = $bp + $np;
	{
	    no warnings qw/uninitialized/;
	    $mass->[$i] = $val;
	}
	if ($i != $current) {
	    $sum += $val;
	}
    }
    my $logQB = log($mass->[$current] / $sum);

    # nCk is proportional to the inverse of k!(n-k)!
    # k! == gamma(k+1)
    my $cFInvLog = lgam($proposal + 1) + lgam($M - $proposal + 1);
    my $cBInvLog = lgam($current + 1) + lgam($M - $current + 1);
    $logQB -= $cFInvLog;
    $logQF -= $cBInvLog;
    return ($proposal, $logQF, $logQB);
}


# slice_sampler1d() implements the univariate "range doubling" slice sampler
# described in Neal (2003) "Slice Sampling", The Annals of Statistics 31(3), 705-767.
# slice_sampler1d() implements a 1-d MCMC slice sampler.
sub sliceSampler1d {
    my ($logF, $x, $minX, $maxX, $w, $nsamples, $maxNFEval) = @_;
    $minX = -9**9**9 unless (defined($minX));
    $maxX = 9**9**9 unless (defined($maxX));
    $w = 0.0 unless (defined($w));
    $nsamples = 1 unless (defined($nsamples));
    $maxNFEval = 20 unless (defined($maxNFEval));

    if ($w <= 0.0) {
	if ($minX > -9**9**9 && $maxX < 9**9**9) {
	    $w = ($maxX - $minX) / 4;
	} else {
	    $w = max((($x < 0.0) ? -$x : $x) / 4, 0.1);
	}
    }

    my $logFx = &$logF($x);
    for (my $sample = 0; $sample < $nsamples; $sample++) {
	my $logY = $logFx + log(gsl_rng_uniform($rng->raw) + 1e-100);

	my $xl = $x - $w * gsl_rng_uniform($rng->raw); # lower bound on slice interval
	my $logFxl = &$logF($xl);
	my $xr = $xl + $w;                             # upper bound on slice interval
	my $logFxr = &$logF($xr);

	while ($logY < $logFxl || $logY < $logFxr) { # doubling procedure
	    if (gsl_rng_uniform($rng->raw) < 0.5) {
		$logFxl = &$logF($xl -= $xr - $xl);
	    } else {
		$logFxr = &$logF($xr += $xr - $xl);
	    }
	}

	my $xl1 = $xl;
	my $xr1 = $xr;
      outer:
	while (1) {                              # shrinking procedure
	    my $x1 = $xl1 + gsl_rng_uniform($rng->raw) * ($xr1 - $xl1);
	    my $flag;
	    if ($logY < &$logF($x1)) {
		my $xl2 = $xl;                   # acceptance procedure
		my $xr2 = $xr;
		my $d = 0;

		while ($xr2 - $xl2 > 1.1 * $w) {
		    my $xm = ($xl2 + $xr2) / 2;
		    if (($x < $xm && $x1 >= $xm) || ($x >= $xm && $x1 < $xm)) {
			$d = 1;
		    }
		    if ($x1 < $xm) {
			$xr2 = $xm;
		    } else {
			$xl2 = $xm;
		    }
		    if ($d && $logY >= &$logF($xl2) && $logY >= &$logF($xr2)) {
			$flag = 'unacceptance';
			last;
			# goto unacceptable;
		    }
		}
		if (!$flag || $flag ne 'unacceptance') {
		    $x = $x1;
		    last outer;
		    # $flag = 'acceptance';
		    # goto acceptable;
		}
	    }
	    if (!$flag || $flag ne 'unacceptance') {
		last outer;
		# goto acceptable;
	    }

	  # unacceptable:
	    if ($x1 < $x) {              # rest of shrinking procedure
		$xl1 = $x1;
	    } else {
		$xr1 = $x1;
	    }
	}
	# acceptable:
	$w = (4 * $w + ($xr1 - $xl1)) / 5;   # update width estimate
    }
    return $x;
}

1;
