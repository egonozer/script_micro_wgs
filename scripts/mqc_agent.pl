#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Spec;

my $dir = File::Spec->rel2abs($ARGV[0]);
my @files = glob($dir . "/*.agent.statistics.txt");

print qq{#section_name: "AGEnt"
#description: "Accessory genome characteristics"
#plot_type: "table"
#pconfig:
#  namespace: "AGEnt"
#headers:
#  total_bp:
#    title: "Accessory size"
#    description: "Total accessory genome size, in bp"
#    suffix: " bp"
#    format: "{:,.0f}"
#  num_segs:
#    title: "AGE #"
#    description: "Number of separate accessory genome elements"
#    format: "{:,.0f}"
#  gc_%:
#    title: "% GC"
#    description: "Percent GC of the accessory genome"
#    suffix: "%"
#    format: "{:.2f}"
#  num_cds:
#    title: "CDS #"
#    description: "Number of coding sequences in the accessory genome"
#    format: "{:,.0f}"
Sample  total_bp    gc_%    num_segs	num_cds
};

my @results;
foreach my $file (@files){
    my $id = basename($file);
    $id =~ s/.agent.statistics.txt//;
    open (my $in, "<$file") or die "ERROR: Can't open $file:$\n";
    while (my $line = <$in>){
        chomp $line;
        if ($line =~ m/^accessory/){
            my @tmp = split("\t", $line);
            print "$id\t",join("\t", @tmp[1,2,3,8]), "\n"; 
        }   
    }
    close ($in);
}
