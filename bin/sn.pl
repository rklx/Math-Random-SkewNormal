#!/usr/bin/env perl6
# For Emacs: -*- mode:cperl; mode:folding; -*-
#
# (c) Jiri Vaclavik
#

use Math::Random::SkewNormal;

sub MAIN(Rat $skewness, Int $number_of_realizations) {
    for 1 .. $number_of_realizations {
        printf("%.15f\n", generate_sn($skewness));
    }
}
