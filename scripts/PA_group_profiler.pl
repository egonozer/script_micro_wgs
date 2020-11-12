#!/usr/bin/perl

use strict;
use warnings;
use List::Util::XS;
use List::Util qw/sum/;

## get the full path to the script directory
use File::Basename;
use Cwd 'abs_path';
my $scriptdir = abs_path(dirname(__FILE__));

my $kfile = "$scriptdir/kmers_and_alleles.k21.txt";

die "ERROR: Can't locate $kfile. Please make sure the file is in the same directory as $0\n" unless -e $kfile;

my $usage = "

PA_group_profiler.pl <file.fasta(.gz)>

Searches for group-defining variants and assigns putative group (A, B, or other) to
a genome sequence. Fasta files can be gzipped.

kmer file:
$kfile


";

die $usage unless @ARGV;
my @files = @ARGV;

## read in kmers
my @klist;
my $arm = 0;
my @kall;
my %kalleles;
my $ptot;
open (my $kin, "$kfile") or die "ERROR: Can't open $kfile: $!\n";
while (my $line = <$kin>){
    chomp $line;
    next if $line =~ m/^#/;
    $ptot++;
    my ($pos, $kmer, $abase, $bbase) = split("\t", $line);
    unless ($arm){
        $arm = (length($kmer) - 1 ) / 2;
    }
    push @kall, $kmer;
    die "ERROR: kmer $kmer at position $pos is repeated\n" if $kalleles{$kmer}; 
    @{$kalleles{$kmer}} = ($abase, $bbase);
    my $rkmer = reverse($kmer);
    $rkmer =~ tr/ACGT/TGCA/;
    my ($rabase, $rbbase) = ($abase, $bbase);
    push @kall, $rkmer;
    die "ERROR: reversed kmer $rkmer at position $pos is repeated\n" if $kalleles{$rkmer}; 
    #push @klist, ([$pos, $kmer, $abase, $bbase, $rkmer]);
}
close ($kin);

my $kstring = join("|", @kall);

print "ID\tGroup\t#A_sites\t#B_sites\t#mixed_sites\t#missing_sites\n";
foreach my $file (@files){
    my $in;
    if ($file =~ m/\.gz$/){
        open ($in, "gzip -cd $file | ");
    } else {
        open ($in, "<$file") or die "ERROR: Can't open $file: $!\n";
    }
    $file =~ s/.*\///g;
    if ($file =~ /\./){
        $file =~ s/\.gz$//;
        $file =~ s/.f[^.]*$//;
    }
    my $seq;
    while (my $line = <$in>){
        if ($line =~ m/^>/){
            if ($seq){
                $seq .= "N";
            }
            next;
        }
        $line =~ s/\s//g;
        $line = uc($line);
        $seq .= $line;
    }
    close ($in);
    #print STDERR "Checking $file for kmers\n";
    (my @hits) = $seq =~ m/(?=($kstring))/g;
    #print STDERR "Screening kmers...\n";
    my %counts;
    foreach my $kmer (@hits){
        my $base = substr($kmer, $arm, 1, ".");
        next if $base eq "N";
        if (!$kalleles{$kmer}){
            $kmer = reverse($kmer);
            $kmer =~ tr/ACGT/TGCA/;
            $base =~ tr/ACGT/TGCA/;
            die "ERROR: kmer $kmer not found in %kalleles\n" unless $kalleles{$kmer};
        }
        unless ($counts{$kmer}){
            @{$counts{$kmer}} = (0) x 3;
        }
        my ($abase, $bbase) = @{$kalleles{$kmer}};
        if ($base eq $abase){
            $counts{$kmer}[0]++;
        } elsif ($base eq $bbase) {
            $counts{$kmer}[1]++;
        } else {
            $counts{$kmer}[2]++;
        }
    }
    my ($acount, $bcount, $ocount, $mixed, $missing) = (0) x 5;
    foreach my $kmer (keys %kalleles){
        if ($counts{$kmer}){
            my @tmp = @{$counts{$kmer}};
            my $sum = sum(@tmp);
            if ($sum == 0){
                $missing++;
                next;
            } elsif ($tmp[0] == $sum){
                $acount++;
                next;
            } elsif ($tmp[1] == $sum){
                $bcount++;
                next;
            } elsif ($tmp[2] == $sum){
                $ocount++;
                next;
            } else {
                $mixed++;
            }
        } else {
            $missing++;
        }
    }
    
    #print STDERR "\n$file,$acount,$bcount,$mixed,$missing\n";
    
    my $group = "A";
    $group = "B" if $bcount > $acount;
    
    my $mfrac = ($mixed + $missing) / $ptot;
    my $smaller = (sort{$a <=> $b}($acount, $bcount))[0];
    my $dfrac = $smaller / ($acount + $bcount);
    ## Built these heuristics by examining patterns 739 NCBI isolates
    if ($mfrac > 0.4 or ($mfrac >= 0.1 and $dfrac > 0.1) or ($mfrac < 0.1 and $dfrac > 30)){
        $group = "?";
    }
    
    
    #if ($mfrac <= 0.05){ #can be missing or mixed at up to 5% of sites
    #    if ($acount / $ptot >= 0.8){
    #        $group = "A";
    #    }
    #    if ($bcount / $ptot >= 0.8){
    #        $group = "B";
    #    }
    #}
    print "$file\t$group\t$acount\t$bcount\t$mixed\t$missing\n";
    
    
}