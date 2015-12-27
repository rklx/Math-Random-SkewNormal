unit module Math::Random::SkewNormal;

my Num $pi = atan2(0, -1);

sub generate_sn(Rat $skewness) is export {
    my Num $a = generate_n();
    my Num $b = generate_n();
    my Num $sconst = $skewness / sqrt(1 + $skewness ** 2);
    my Num $c = $sconst * $a + sqrt(1 - $sconst ** 2) * $b;
    return $a > 0 ?? $c !! -$c;
}

sub generate_sn_multi(List $delta, List $omega) is export {
    my $n = @$delta;
    return if $n != @$omega;

    # build correlation matrix
    my $cor = [];
    for 0 .. $n - 1 {
        $cor.[$_][$_] = 1;
    }
    for 1 .. $n {
        $cor.[0][$_] = $delta.[$_-1];
        $cor.[$_][0] = $delta.[$_-1];
    }
    for 1 .. $n -> $i {
        for 1 .. $n -> $j {
            $cor.[$i][$j] = $omega.[$i-1][$j-1];
        }
    }

    my $sample = generate_n_multi($cor);

    if (@$sample.shift > 0) {
       @$sample = map {-$_}, @$sample;
    }

    return $sample;
}

sub generate_n {
    my $u = generate_r();
    my $v = generate_r();
    return sqrt(-2 * log $u) * cos(2 * $pi * $v);
}

sub generate_r {
    return rand;
}

sub generate_n_multi(List $cov) {
    my $n = @$cov || return;
    my $independent = [];
    my $result = [];
    my $sqrt = _choleski($cov) // return;
    for 0 .. $n - 1 {
        $independent.[$_] = generate_n();
    }
    for 0 .. $n - 1 -> $i {
        $result.[$i] = 0;
        for 0 .. $n - 1 -> $j {
            $result.[$i] += $independent.[$j] * $sqrt.[$i][$j];
        }
    }

    return $result;
}

# choleski's decomposition
sub _choleski (List $A) {
    my $L = [];

    for 0 .. @$A - 1 -> $c {
        my $sum = 0;
        for 0 .. $c - 1 -> $j {
            my $i = $c - 1 - $j;
            $sum += $L.[$c][$i] ** 2;
        }
        $L.[$c][$c] = sqrt($A.[$c][$c] - $sum);
        for $c + 1 .. @$A - 1 -> $r {
            my $sum = 0;
            for 0 .. $c - 1 -> $j {
                my $i = $c - 1 - $j;
                $sum += $L.[$r][$i] * $L.[$c][$i];
            }
            $L.[$r][$c] = ($A.[$r][$c] - $sum) / $L.[$c][$c];
        }
        for 0 .. $c - 1 -> $r {
            $L.[$r][$c] = 0;
        }
    }

    return $L;
}


=begin pod

=head1 NAME

Math::Random::SkewNormal - Handy, easy-to-use Skew-Normal random number generator

=head1 SYNOPSIS

  use Math::Random::SkewNormal;
  
  # realizations for skew normal distribution (special case):
  $skewness = 1;
  generate_sn($skewness);

  # realizations for multidimensional skew normal distribution:
  $delta = [0, 0];
  $omega = [[1, 0], [0, 1]];
  generate_sn_multi($delta, $omega);

=head1 DESCRIPTION

This module transforms uniformly spaced random variable realizations into
realizations that follow the Skew-Normal (I<SN>) Probability Density Function (I<PDF>).

We accept following definition of the Skew-Normal Distribution:

=over

=item 1-dimensional SN is determined by PDF: f(x,a) = 2 phi(x) Phi(a x),
    where phi is the PDF and Phi is the CDF of Normal distribution

=item Let X = (X_0, ..., X_k) follows distribution N_k+1 (0, Omega*), where Omega* is
    the matrix of the form:

       [   1   delta ]
       [ delta Omega ]

Omega is a matrix and delta is a vector. Then (X_1, ..., X_n)|(X_0 > 0) follows
n-dimensional Skew Normal distribution.


=back

=head1 ALGORITHM

We use following algorithms:

=over

=item * Box-Muller transformation; for generation Normal realizations

=item * Choleski's decomposition; for inverting matrices

=item * A lemma mentioned e. g. in article from A. Azzalini:
    The skew-normal distribution and related multivariate families;
    for generating Skew-Normal realizations

=back

=head1 FUNCTIONS

=head2 generate_sn

Returns realization of SN(0,1, a), where B<a> is the only parameter.
The meaning of B<a> is "skewness" (see definition above).

=head2 generate_sn_multi

Returns realization of centralized n-dimensional Skew-Normal distribution.

1st parameter is correlation matrix of the form C<[[row_1], [row_2], ..., [row_n]]>,
for example C<[[1, 0], [0, 1]]> for 2-dimensional unit matrix.

2nd parameter is delta vector of the form C<[delta_1, ..., delta_n]>, e. g. C<[0, 1]>.

=head1 SEE ALSO

L<Math::Random>, L<Math::Random::Cauchy>,

The examples in F<examples/> subdirectory of this distribution
contains example of generated data with histogram.

The scripts in F<bin/> directory of this distribution contains
tools for generating realizations on command line.

=head1 AUTHOR

Jiri Vaclavik, E<lt>my name dot my last name at gmail dot comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Jiri Vaclavik

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl 6 itself.

=end pod
