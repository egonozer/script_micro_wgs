#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
our ($opt_d);
getopts('d:');

my $dir = $opt_d;

foreach my $sample (@ARGV){
    print "$dir/results/agent/$sample.agent.accessory.fasta\t$sample\n";
    print STDERR "$dir/results/agent/$sample.agent.accessory_loci.txt\t$sample\n";
}
