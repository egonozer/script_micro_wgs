#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my $usage = "
snippy_process.pl

Takes output from Snippy and produces two files:
    1) snps.processed_nt_variants.tsv:
        All nucleotide SNP variants
    2) snps.processed_aa_variants.tsv:
        Amino acid changes for all non-synonymous variants
    Both will be output in the same directory as the input snps.tab file
        
Required:
  -s    snps.tab file output by Snippy
  
Optional:
  -b    bed file of targets to restrict output to (i.e. core genome)
  -d    debug mode

";



use Getopt::Std;
our ($opt_s, $opt_b, $opt_d);
getopts('s:b:d');

die $usage unless $opt_s;

my $file    = $opt_s;
my $bed     = $opt_b if $opt_b;

$file = abs_path($file);
my ($name, $path) = fileparse($file);

my %bedhash;
if ($bed){
    open (my $in, "<$bed") or die "ERROR: Can't open $bed: $!\n";
    while (my $line = <$in>){
        chomp $line;
        my ($chrom, $start, $stop) = split("\t", $line);
        $start++;
        for my $i ($start .. $stop){
            $bedhash{$chrom}{$i} = 1;
        }
    }
    close ($in);
}

open (my $in, "<$file") or die "ERROR: Can't open $file:\n";
my ($nfile, $afile) = ("$path/snps.processed_nt_variants.tsv", "$path/snps.processed_aa_variants.tsv");
($nfile, $afile) = ("snps.processed_nt_variants.tsv", "snps.processed_aa_variants.tsv") if $opt_d;
open (my $nout, ">$nfile");
print $nout "CHROM\tPOS\tREF\tALT\n";
open (my $aout, ">$afile");
print $aout "CHROM\tLOCUS_ID\tPOS\tREF\tALT\n";
while (my $line = <$in>){
    chomp $line;
    my @tmp = split("\t", $line);
    my ($chrom, $pos, $type, $ref, $alt, $ftype, $effect, $loc, $gene, $prod) = @tmp[0,1,2,3,4,6,10,11,12,13];
    if ($bed){
        next unless $bedhash{$chrom}{$pos};
    }
    if ($type =~ m/snp|complex|mnp/){
        ## for unequal length ref/alt that are complex variants that include indels, we'll just check the first base and count it as a SNP if it's different. 
        ## Everyting else in the uneven complex variants is just too difficult to reconcile.
        if (length($ref) != length($alt)){
            $ref = substr($ref, 0, 1);
            $alt = substr($alt, 0, 1);
        }
        ## get the variant positions
        my $mask = $ref ^ $alt;
        while ($mask =~ m/[^\0]/g){
            my $mpos = $-[0];
            my $refbase = substr($ref, $mpos, 1);
            my $altbase = substr($alt, $mpos, 1);
            my $opos = $pos + $mpos;
            print $nout "$chrom\t$opos\t$refbase\t$altbase\n";
        }
    }
    ## evaluate amino acid changes
    ## "X" = frameshift or start codon loss
    next unless $effect;
    (my $var) = $effect =~ m/p\.(\D+\d+\D+)/;
    next unless $var;
    my ($ref_a, $pos_a, $alt_a) = split(/(\d+)/, $var) if $var;
    next unless ($ref_a ne $alt_a);
    ## Will ignore insertions for this analysis
    if (substr($alt_a,0,1) ne uc(substr($alt_a,0,1))){
        if ($alt_a eq "fs"){ #frameshift
            my $out_r = oneletter($ref_a);
            print $aout "$chrom\t$loc\t$pos_a\t$out_r\tX\n";
        }
        ## otherwise will skip "ins", "del", and "dup"
    } else {
        if ($alt_a eq "?"){ ## start codon loss
            my $out_r = oneletter($ref_a);
            print $aout "$chrom\t$loc\t$pos_a\t$out_r\tX\n";
        }
        if ($alt_a eq "*") { ## premature stop codon
            my $out_r = oneletter(substr($ref_a,0,3));
            print $aout "$chrom\t$loc\t$pos_a\t$out_r\t*\n";
        }
        if ($alt_a =~ m/ext\*\?$/){ ## stop lost
            $alt_a =~ s/ext\*\?$//; #remove the 'ext*?'. Next if statement will handle the amino acid changes. 
        }
        if (length($ref_a) == length($alt_a)){
            if (length($ref_a) == 3){
                my ($out_r, $out_a) = (oneletter($ref_a), oneletter($alt_a));
                if ($out_r eq "X" or $out_a eq "X"){
                    die "ERROR: undefined amino acid in snippy output: $line\n"; ## Haven't seen this yet.
                }
                print $aout "$chrom\t$loc\t$pos_a\t$out_r\t$out_a\n";
            } else {
                my @rarray = ($ref_a =~ m/(.{3})/g);
                my @aarray = ($alt_a =~ m/(.{3})/g);
                for my $i (0 .. $#rarray){
                    next if $rarray[$i] eq $aarray[$i];
                    my $opos = $pos_a + $i;
                    my ($out_r, $out_a) = (oneletter($rarray[$i]), oneletter($aarray[$i]));
                    if ($out_r eq "X" or $out_a eq "X"){
                        die "ERROR: undefined amino acid in snippy output: $line\n"; ## Haven't seen this yet.
                    }
                    print $aout "$chrom\t$loc\t$opos\t$out_r\t$out_a\n";
                }
            }
        }
    }
}

sub oneletter {
    my $input = shift;
    my %hash = (
        'Ala'=>'A',
        'Arg'=>'R',
        'Asn'=>'N',
        'Asp'=>'D',
        'Cys'=>'C',
        'Gln'=>'Q',
        'Glu'=>'E',
        'Gly'=>'G',
        'His'=>'H',
        'Ile'=>'I',
        'Leu'=>'L',
        'Lys'=>'K',
        'Met'=>'M',
        'Phe'=>'F',
        'Pro'=>'P',
        'Ser'=>'S',
        'Thr'=>'T',
        'Trp'=>'W',
        'Tyr'=>'Y',
        'Val'=>'V',
        'Ter'=>'*'
    );
    my $aa = "X";
    $aa = $hash{$input} if $hash{$input};
    return($aa);
}
