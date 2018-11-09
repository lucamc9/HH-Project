#!/usr/bin/perl
use warnings;
use strict;
use integer;
my %a;
while (<>) {
   my ($gene, $n) = /^\s*(\S+)\s*\t\s*(\S+)/;
      $a{$gene} += $n if defined $n;
}
print "$_\t${a{$_}}\n" for sort keys %a;
