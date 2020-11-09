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

";

use Getopt::Std;
our ($opt_s, $opt_b);
getopts('s:b:');

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
open (my $nout, ">$path/snps.processed_nt_variants.tsv");
print $nout "CHROM\tPOS\tREF\tALT\n";
open (my $aout, ">$path/snps.processed_aa_variants.tsv");
print $aout "CHROM\tLOCUS_ID\tPOS\tREF\tALT\n";
while (my $line = <$in>){
    chomp $line;
    my @tmp = split("\t", $line);
    my ($chrom, $pos, $type, $ref, $alt, $ftype, $effect, $loc) = @tmp[0,1,2,3,4,6,10,11];
    if ($bed){
        next unless $bedhash{$chrom}{$pos};
    }
    next unless $type =~ m/snp|complex/;
    if ($type eq "snp"){
        print $nout "$chrom\t$pos\t$ref\t$alt\n";
    } else { ## complex
        ## trim down unequal length sequences before pulling out variants
        my $maxlength = (sort{$a <=> $b}(length($ref), length($alt)))[0];
        $ref = substr($ref, 0, $maxlength);
        $alt = substr($alt, 0, $maxlength);
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
    next unless $effect;
    (my $var) = $effect =~ m/p\.(\D+\d+\D+)/;
    my ($ref_a, $pos_a, $alt_a) = split(/(\d+)/, $var) if $var;
    if ($effect =~ m/missense_variant/){
        if (length($ref_a) == length($alt_a)){
            if (length($ref_a) == 3){
                my ($out_r, $out_a) = (oneletter($ref_a), oneletter($alt_a));
                print $aout "$chrom\t$loc\t$pos_a\t$out_r\t$out_a\n";
            } else {
                my @rarray = ($ref_a =~ m/(.{3})/g);
                my @aarray = ($alt_a =~ m/(.{3})/g);
                for my $i (0 .. $#rarray){
                    next if $rarray[$i] eq $aarray[$i];
                    my $opos = $pos_a + $i;
                    next unless $bedhash{$chrom}{$opos};
                    my ($out_r, $out_a) = (oneletter($rarray[$i]), oneletter($aarray[$i]));
                    print $aout "$chrom\t$loc\t$opos\t$out_r\t$out_a\n";
                }
            }
        } else {
            ### I don't know how best to deal with uneven complex variants and they seem relatively rare, so will ignore for now. Sorry. 
        }
    }
    
    next unless $type eq "snp"; ## for now. Need some examples to deal with stop codons in complex
    if ($effect =~ m/stop_gained/){
        (my $var) = $effect =~ m/p.(\D+\d+\*)/;
        my ($ref_a, $pos_a, $alt_a) = split(/(\d+)/, $var) if $var;
        my $out_r = oneletter($ref_a);
        print $aout "$chrom\t$loc\t$pos_a\t$out_r\t$alt_a\n";
    }
    if ($effect =~ m/stop_lost/){
        $alt_a = substr($alt_a, 0, 3); ## stop_lost always annotated with tailing "ext*?" for splice_region_variant
        my ($out_r, $out_a) = (oneletter($ref_a), oneletter($alt_a));
        print $aout "$chrom\t$loc\t$pos_a\t$out_r\t$out_a\n";
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
