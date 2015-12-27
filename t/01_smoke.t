#!/usr/bin/env perl6
use v6;
use Test;
use lib 'lib';

use Math::Random::SkewNormal;

my Rat $skewness = 1.0;
my List $delta    = [0, 0];
my List $omega    = [[1, 0], [0, 1]];

# following functions generate realizations:
is generate_sn($skewness).WHAT, Num,  'generate_sn: is a number';
is generate_sn_multi($delta, $omega).end + 1, 2, 'generate_sn_multi: got vector of correct size';

done-testing();