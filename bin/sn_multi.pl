#!/usr/bin/env perl
# For Emacs: -*- mode:cperl; mode:folding; -*-
#
# (c) Jiri Vaclavik
#
use Math::Random::SkewNormal;

sub MAIN(Int :$number_of_realizations, List :$delta, List :$omega_matrix_by_rows) {
	my Int $dim = $delta.end + 1;
	my List $A;

    for 0 .. $dim - 1 -> $row {
        for 0 .. $dim - 1 -> $col {
            $A.[$row][$col] = $omega_matrix_by_rows.shift;
        }
    }

    for 1 .. $number_of_realizations {
        print @(generate_sn_multi($delta, $A)), "\n";
    }
}
