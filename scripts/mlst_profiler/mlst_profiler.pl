#!/usr/bin/perl

my $version = "0.3.1";

## Changes from v0.3
## Make allele sequence files more portable by allowing relative paths
## Use blast in PATH or designate another blast
## Remove makeblastdb dependency. Becomes problematic when running on a cluster.

## Changes from v0.2.1
## Blast against all MLST targets at once, rather than one at a time. Should be MUCH faster.

## Changes from v0.2
## Modified to be more flexible and be able to do profiling from other PubMLST allele/profile combinations

## Changes from v0.1
## Automatically correct for issue in PA acs allele definitions where alleles 38 and 172 are both present in isolates with ST235.
## Make blast location modifiable from command line
## Make alignment slop user-definable

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions qw ( catfile path );


$|++;

my $blastdir = "NA";
my $path = is_path("blastn");
if ($path){
    $blastdir = $path;
}

my $usage = "
mlst_profiler.pl (version $version)

Determines the MLST profile from a fasta sequence.

Required:
  -p    File with MLST profiles. Format should be the same as profile files
        downloaded from pubmlst.org, i.e. a header line with the allele names
        starting with the sequence type ID, separated by tabs. Following lines
        list sequence types and allele IDs, separated by tabs.
        
        Example:
        ST	adk	atpA	dxr	glyA	recA	sodA	tpi
        1	1	1	1	10	1	3	5
        2	1	1	2	1	5	3	1
        3	1	1	2	1	1	1	1
        ...
        
  -a    File with list of allele sequence files, one per line.
        Format should be: allele name (matching the column headers in the MLST
        profile above), tab, path to fasta file with allele sequences. Allele
        headers should be the allele name, underscore, allele number
        (for example: adk_1 is the first variant of the adk allele).
        File example:
        adk /path/to/adk.fas
        atpA    /path/to/atpA.fas
        ...

  -f    Query fasta sequence file. Can be gzipped.
    OR
  -F    list of query fasta sequence files, one path to sequence file per line.
        Can be gzipped.
        Optional: sequence ID can follow the path, separated by a tab.
          Example: /path/to/file.fasta<tab>genome1

Options:
  -s    slop, i.e. number of bases by which a potential novel allele can differ
        from its most closely related allele and not be considered potentially
        truncated
        (default: 2)
  -n    prefix for sequences of potential novel alleles. If this is given, will
        keep track of potential novel alleles not found in the given allele
        sequence files and output their sequences.
        (default: novel allele sequences will not be tracked)

  -b    path to Blast+ binaries
        (default: '$blastdir')

";

use Getopt::Std;
use vars qw( $opt_p $opt_a $opt_f $opt_F $opt_b $opt_s $opt_n );
getopts('p:a:f:F:b:s:n:');

die $usage unless ($opt_p and $opt_a and ($opt_f or $opt_F));
die $usage if ($opt_f and $opt_F);
die $usage if ($blastdir eq "NA" and !$opt_b);
my $pfile       = $opt_p;
my $afile       = $opt_a;
my $qfile       = $opt_f if $opt_f;
my $qlistfile   = $opt_F if $opt_F;
my $slopleng    = $opt_s ? $opt_s : 2;
$blastdir       = $opt_b ? $opt_b : $blastdir;
my $npref       = $opt_n if $opt_n;

# Read in the allele file list and create a temporary sequence file with all alleles
$afile = abs_path($afile);
my ($aname, $apath) = fileparse($afile);
open (my $ain, "<$afile") or die "ERROR: Can't open $afile: $!\n";
#open (my $aout, ">temp_alleles.fasta");
my %afiles;
my @apaths;
my $a_tot_count = 0;
while (my $line = <$ain>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    $line =~ s/\s*$//g;
    my ($fid, $fpath) = split("\t", $line);
    $fid =~ s/\s//g;
    
    ## Check to see if the filepath is absolute or relative:
    unless (-e $fpath){
        my $temppath = abs_path("$apath/$fpath");
        if (-e "$temppath"){
            $fpath = $temppath;
        } else {
            die "ERROR: Can't locate sequence file $fpath for allele $fid\n";
        }
    }
    
    ## add all paths to an array rather than readig them into a temporary file
    push @apaths, $fpath;
    
    open (my $tin, "<$fpath") or die "ERROR: Can't open $path: $!\n";
    while (my $tline = <$tin>){
        #print $aout $tline;
        $a_tot_count++ if $tline =~ m/^>/;
    }
    close($tin);
    $afiles{$fid} = 1;
}
close ($ain);
my $apath_string = join(" ", @apaths);
#close ($aout);
#my $status = system("$blastdir/makeblastdb -in temp_alleles.fasta -dbtype nucl >/dev/null");
#die "ERROR: makeblastdb exited with status $status\n" if $status;

## multiply the a_tot_count by 100, just to try to be sure to get all alignments
$a_tot_count = $a_tot_count * 100;

## read the profile matrix into a hash
#open (my $pin, "<$pfile") or die "ERROR: Can't open $pfile: $!\n";
#my %profile;
#my @allele_order;
#my %seen_alleles;
#my %mlst_clade; #for C. diff
#my $mc; #for C. diff
#while (my $line = <$pin>){
#    chomp $line;
#    next if $line =~ m/^\s*$/;
#    my @tmp = split("\t", $line);
#    my $st = shift @tmp;
#    if (!@allele_order){
#        for my $i (0 .. $#tmp){
#            my $allele = $tmp[$i];
#            unless ($afiles{$allele}){
#                $mc = $i if $allele eq "mlst_clade";
#                next;
#            }
#            $seen_alleles{$allele}++;
#            push @allele_order, $allele;
#        }
#        next;
#    }
#    my @vals;
#    for my $i (0 .. $#allele_order){
#        push @vals, $tmp[$i];
#    }
#    my $val_string = join(",", @vals);
#    die "ERROR: profile $val_string already exists as ST $profile{$val_string}\n" if $profile{$val_string};
#    $profile{$val_string} = $st;
#    $mlst_clade{$val_string} = $tmp[$mc] if defined $mc and $tmp[$mc];
#}
#close ($pin);
#foreach my $al (sort keys %afiles){
#    die "ERROR: Allele definition $al was not found in the profile file \"$pfile\"\n" unless $seen_alleles{$al};
#}

# read the alternate profile matrix into a hash
my %profile;
my @allele_order;
my %seen_alleles;
my $header;
open (my $pin, "<$pfile") or die "ERROR: Can't open $pfile: $!\n";
my %non_alleles;
while (my $line = <$pin>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    my @tmp = split("\t", $line);
    if (!@allele_order){
        for my $i (0 .. $#tmp){
            my $allele = $tmp[$i];
            unless ($afiles{$allele}){
                $non_alleles{$i} = 1;
                $header .= "$allele\t";
                next;
            }
            $seen_alleles{$allele}++;
            push @allele_order, $allele;
        }
        chop $header; #remove final tab character
        next;
    }
    my @vals;
    my $st;
    for my $i (0 .. $#tmp){
        if ($non_alleles{$i}){
            $st .= "$tmp[$i]\t";
            next;
        }
        push @vals, $tmp[$i];
    }
    chop $st; #remove final tab character;
    my $val_string = join(",", @vals);
    die "ERROR: profile $val_string (ST $st) already exists as ST $profile{$val_string}\n" if $profile{$val_string};
    $profile{$val_string} = $st;
}
close ($pin);
foreach my $al (sort keys %afiles){
    die "ERROR: Allele definition $al was not found in the profile file \"$pfile\"\n" unless $seen_alleles{$al};
}
my $non_count = scalar keys %non_alleles;

print "file\t", join("\t", @allele_order), "\t$header\n"; #header
#print "file\t", join("\t", @allele_order), "\tST"; #header
#print "\tclade" if %mlst_clade;
#print "\n";

my @qlist;
if ($opt_F){
    open (my $fin, "<$qlistfile") or die "ERROR: Can't open $qlistfile: $!\n";
    while (my $line = <$fin>){
        chomp $line;
        next if $line =~ m/^\s*$/;
        my ($path, $id) = split("\t", $line);
        $path =~ s/\s*$//;
        $path =~ s/^\s*//;
        unless ($id){
            $id = $path;
            $id =~ s/\/*[^\/]+\///g; #remove path
            $id =~ s/\.gz$//; #remove gz
            $id =~ s/\.[^.]+$//; #remove last suffix
        }
        $path =~ s/\|/\\\|/;
        push @qlist, ([$path, $id]);
    }
} else {
    my $id = $qfile;
    $id =~ s/\/*[^\/]+\///g; #remove path
    $id =~ s/\.gz$//; #remove gz
    $id =~ s/\.[^.]+$//; #remove last suffix
    $qfile =~ s/\|/\\\|/;
    push @qlist, ([$qfile, $id]);
}

my %novel_alleles;
my %nfiles;
foreach my $slice (@qlist){
    my ($file, $fileid) = @{$slice};
    print STDERR "PROFILING $fileid...\n";
    
    # If the query file is gzipped, unzip it to a temporary file
    my $isgzip;
    if ($file =~ m/\.gz$/){
        system("gzip -cd $file > qtemp.fasta");
        $file = "qtemp.fasta";
        $isgzip = 1;
    }
    
    # Blast against the alleles to get the allele IDs
    print STDERR "Blasting...\n";
    my %hits;
    open (my $bin, "cat $apath_string | $blastdir/blastn -query $file -subject - -max_target_seqs $a_tot_count -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send pident' | ") or die "ERROR: blastn failed: $!\n";
    
    #open (my $bin, "$blastdir/blastn -query $file -db temp_alleles.fasta -max_target_seqs $a_tot_count -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send pident' | ") or die "ERROR: blastn failed: $!\n";
    while (my $line = <$bin>){
        chomp $line;
        my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $pident) = split("\t", $line);
        (my $aid) = $sseqid =~ m/^(\S+)_\d+$/;
        push @{$hits{$aid}}, $line;
    }
    close ($bin);
    
    # Determine alleles.
    my @allele_vals;
    foreach my $allele (@allele_order){
        die "ERROR: No allele sequence file given for $allele\n" unless $afiles{$allele};
        print STDERR "Determining $allele...\n";
        my $allele_num;
        my ($top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line);
        my $trunc_hit;
        my $out_allele = "-";
        if ($hits{$allele}){
            foreach my $line (@{$hits{$allele}}){
                my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $pident) = split("\t", $line);
                my $aleng = abs($send - $sstart) + 1;
                (my $anum) = $sseqid =~ m/_(\d+)$/;
                if ($pident == 100){
                    if ($aleng == $slen){
                        if (!$allele_num){
                            ($allele_num, $top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line) = ($anum, $anum, $qseqid, $qstart, $qend, $pident, $line);
                            next;
                        } else {
                            if ($anum == $allele_num){
                                print STDERR "Two identical alleles for $allele exist:\n\t$line\n";
                                next;
                            } else {
                                if ($allele_num =~ m/N/){
                                    ($allele_num, $top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line) = ($anum, $anum, $qseqid, $qstart, $qend, $pident, $line);
                                } else {
                                    print STDERR "At least two different alleles for $allele exist:\n\t$line\n";
                                    my @tmp = split("&", $allele_num);                          
                                    push @tmp, $anum;
                                    @tmp = grep {s/(^|\D)0+(\d)/$1$2/g,1} sort grep {s/(\d+)/sprintf"%09.9d",$1/ge,1} @tmp; #natural sort. See www.perlmonks.org/?node_id=68185
                                    $allele_num = join("&", @tmp);
                                }
                                next;
                            }
                        }
                    } else {
                        print STDERR "Allele $allele is truncated\n\t$line\n" unless $allele_num;
                        $trunc_hit = 1;
                    }
                } else {
                    if ($top_qid){
                        next if ($qseqid eq $top_qid and $qstart <= $top_qend and $qend >= $top_qstart); #skip if this overlaps a previous hit
                        #print STDERR "Second sub-ideal hit found: $line\n"; ## Not worth keeping track of if there's already a top hit
                    } else {
                        next if $trunc_hit;
                        if (slop($slen,$aleng,$slopleng)){ #hit can differ by up to 2 bases from allele length;
                            ($top_qid, $top_qstart, $top_qend, $top_line) = ($qseqid, $qstart, $qend, $line);
                            print STDERR "Potential novel allele found: $line\n";
                            
                            # get novel allele sequence
                            if (defined $npref){
                                open (my $nin, "<$file");
                                my $seq;
                                my $readseq;
                                while (my $line = <$nin>){
                                    chomp $line;
                                    next if $line =~ m/^\s.*$/;
                                    if ($line =~ m/^>/){
                                        last if $readseq;
                                        my $id = substr($line, 1);
                                        $id =~ s/\s.*$//;
                                        $readseq = 1 if $id eq $qseqid;
                                        next;
                                    }
                                    next unless $readseq;
                                    $line =~ s/\s//g;
                                    $seq .= $line;
                                }
                                close ($nin);
                                if ($seq){
                                    my $dist = $qend - $qstart + 1;
                                    my $aseq = substr($seq, $qstart - 1, $dist);
                                    if ($sstart > $send){
                                        $aseq = reverse $aseq;
                                        $aseq =~ tr/ACGTacgt/TGCAtgca/;
                                    }
                                    my $num = 0;
                                    $num = $novel_alleles{$allele}{"num"} if $novel_alleles{$allele}{"num"};
                                    if ($novel_alleles{$allele}{$aseq}){
                                        $num = $novel_alleles{$allele}{$aseq};
                                    } else {
                                        $num++;
                                        $novel_alleles{$allele}{$aseq} = $num;
                                        $novel_alleles{$allele}{"num"} = $num;
                                        my $nout;
                                        if ($nfiles{$allele}){
                                            my $nfile = $nfiles{$allele};
                                            open ($nout, ">>$nfile");
                                        } else {
                                            my $nfile = "$npref.$allele.fasta";
                                            open ($nout, ">$nfile");
                                        }
                                        print $nout ">$allele\_N$num\n$aseq\n";
                                        close ($nout);
                                    }
                                    $allele_num = "N$num";
                                } else {
                                    die "ERROR: Problem trying to find sequence of novel allele: $line\n";
                                }
                            }
                            
                        } else {
                            print STDERR "Potential truncated novel allele found: $line\n";
                            $out_allele = "NT";
                        }
                    }
                }
            }
            if ($allele_num){
                $out_allele = $allele_num;
            } else {
                if ($top_qid){
                    $out_allele = "N";
                } elsif (!$top_qid and $trunc_hit) {
                    $out_allele = "T";
                } else {
                    print STDERR "No hits\n";
                }
            }
            push @allele_vals, $out_allele;
            print STDERR "Top line: $top_line\n" if $top_line;
        } else {
            print STDERR "No hits\n";
	    push @allele_vals, $out_allele;
        }
    }
    
    if ($isgzip){
        unlink("$file");
    }
    
    # Determine sequence type
    my @tmpst = ("?") x $non_count;
    my $st = join("\t", @tmpst);
    my $pstring = join(",", @allele_vals);
    if ($profile{$pstring}){
        $st = $profile{$pstring};
    }
    
    # PA ST235 hack
    if ($st eq "?" and ($pstring eq "38&172,11,3,13,1,2,4" or $pstring eq "172&38,11,3,13,1,2,4")){
        $st = "235";
    }
    
    #if (%mlst_clade){
    #    $st .= "\t";
    #    if (my $mc = $mlst_clade{$pstring}){
    #        $st .= "$mc";
    #    }
    #}
    print "$fileid\t", join("\t", @allele_vals), "\t$st\n";
    print STDERR "\n";
}

#foreach my $allele (@allele_order){
#    unlink("tmp_$allele.nhr");
#    unlink("tmp_$allele.nin");
#    unlink("tmp_$allele.nsq");
#}

#my @files = <temp_alleles.fasta*>;
#unlink @files;

#---------------------
sub slop{
    my $val1 = shift;
    my $val2 = shift;
    my $pm = shift;
    return(1) if abs($val1 - $val2) <= $pm;
    return(0);
}

sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($dir) if -x $f;
        }
    }
    return();
}
