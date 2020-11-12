#!/usr/bin/env perl

use strict;
use warnings;
use File::Find;

my $usage = "
snippy_compile.pl -o <output_prefix> </path/to/folder/containing/snippy_output_directories/>

Compiles all of the results from snippy_process.pl into two files:
  <prefix>.processed_aa_variants.tsv
  <prefix>.processed_nt_variants.tsv

";

use Getopt::Std;
use vars qw( $opt_o );
getopts('o:');

die $usage unless @ARGV;

my $pref    = $opt_o ? $opt_o : "snippy";
my $dir     = $ARGV[0];

my @ids;
find(\&dir_search, $dir);
@ids = sort(@ids);

## start with the nt resuts;
my @chrom_order;
my %chrom_seen;
my @nt_ids;
my %hash;
my %nt_ref;
foreach my $id (@ids){
    if (-e "$dir/$id/snps.processed_nt_variants.tsv"){
        push @nt_ids, $id;
        open (my $in, "<$dir/$id/snps.processed_nt_variants.tsv");
        while (my $line = <$in>){
            chomp $line;
            next if $line =~ m/^CHROM/;
            my ($chrom, $pos, $ref, $alt) = split("\t", $line);
            unless ($chrom_seen{$chrom}){
                push @chrom_order, $chrom;
                $chrom_seen{$chrom} = 1;
            }
            $nt_ref{$chrom}{$pos} = $ref;
            $hash{$chrom}{$pos}{$id} = $alt;
        }
        close ($in);
    }
}

if (@nt_ids){
    open (my $out_nt, ">$pref.processed_nt_variants.tsv");
    print $out_nt "CHROM\tPOS\tREF\t";
    print $out_nt join("\t", @nt_ids), "\n";
    foreach my $chrom (@chrom_order){
        foreach my $pos (sort {$a <=> $b} keys %{$hash{$chrom}}){
            my $ref = $nt_ref{$chrom}{$pos};
            print $out_nt "$chrom\t$pos\t$ref";
            foreach my $id (@nt_ids){
                my $val = $ref;
                $val = $hash{$chrom}{$pos}{$id} if $hash{$chrom}{$pos}{$id};
                print $out_nt "\t$val";
            }
            print $out_nt "\n";
        }
    }
    close ($out_nt);
} else {
    print STDERR "WARNING: No processed_nt_variants.tsv files found\n";
}

## now the aa results
%hash = ();
my @aa_ids;
my %aa_ref;
foreach my $id (@ids){
    if (-e "$dir/$id/snps.processed_aa_variants.tsv"){
        push @aa_ids, $id;
        open (my $in, "<$dir/$id/snps.processed_aa_variants.tsv");
        while (my $line = <$in>){
            chomp $line;
            next if $line =~ m/^CHROM/;
            my ($chrom, $lid, $pos, $ref, $alt) = split("\t", $line);
            unless ($chrom_seen{$chrom}){
                push @chrom_order, $chrom;
                $chrom_seen{$chrom} = 1;
            }
            $aa_ref{$chrom}{$lid}{$pos} = $ref;
            $hash{$chrom}{$lid}{$pos}{$id} = $alt;
        }
        close ($in);
    }
}

if (@aa_ids){
    open (my $out_aa, ">$pref.processed_aa_variants.tsv");
    print $out_aa "CHROM\tLOCUS_ID\tPOS\tREF\t";
    print $out_aa join("\t", @aa_ids), "\n";
    foreach my $chrom (@chrom_order){
        foreach my $lid (sort{$a cmp $b} keys %{$hash{$chrom}}){
            foreach my $pos (sort {$a <=> $b} keys %{$hash{$chrom}{$lid}}){
                my $ref = $aa_ref{$chrom}{$lid}{$pos};
                print $out_aa "$chrom\t$lid\t$pos\t$ref";
                foreach my $id (@aa_ids){
                    my $val = $ref;
                    $val = $hash{$chrom}{$lid}{$pos}{$id} if $hash{$chrom}{$lid}{$pos}{$id};
                    print $out_aa "\t$val";
                }
                print $out_aa "\n";
            }
        }
    }
    close ($out_aa);
} else {
    print STDERR "WARNING: No processed_aa_variants.tsv files found\n";
}

#print join("\n", @aa_files), "\n";

sub dir_search {
    my $id = $_;
    push @ids, $id if -e "$dir/$id/snps.processed_nt_variants.tsv";
}

