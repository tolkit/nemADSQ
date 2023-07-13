#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;


# globals

my $EXE = basename($0);
my $VERSION = "0.0.1";
my $AUTHOR = 'Torsten Seemann (@torstenseemann) and Pablo Gonzalez de la Rosa';
my $HOMEPAGE = "undefined";

# SAM file TSV columns
use constant {
  SAM_FLAG  => 1,
  SAM_RNAME => 2,
  SAM_POS   => 3,
  SAM_MAPQ  => 4,
  SAM_CIGAR => 5,
  SAM_TLEN  => 8,
  SAM_SEQ   => 9,
};

#----------------------------------------------------------------------
# command line parameters

my $telomere  = "TTAGGC";
my $ref       = '';
my $ref_index = '';
my $debug     = 0;

#----------------------------------------------------------------------
sub usage {
  my($exitcode) = @_;
  $exitcode=0 if !defined($exitcode) or $exitcode eq 'help';
  my $fh = $exitcode ? \*STDERR : \*STDOUT;
  print $fh
    "SYNOPSIS\n  Add telomeric sequence to reference contigs from clipped alignments\n",
    "AUTHOR\n  $AUTHOR\n",
    "USAGE\n",
    "  % minimap2 -a ref.fa telomeric_reads.fa | teloextend --ref ref.fa > telo_extended.fasta\n",
    "OPTIONS\n",
    "  --help           This help\n",
    "  --version        Print version and exit\n",
    "  --ref FASTA      Reference genome - needs FASTA.fai index\n",
    "  --telomere STR   Telomeric repeat for species (default is nematodes' $telomere)\n",
    "  --debug          Print verbose debug info to stderr\n",
    "HOMEPAGE\n  $HOMEPAGE\n",
    "";
  exit($exitcode);
}

#----------------------------------------------------------------------
# getopts

@ARGV or usage(1);

GetOptions(
  "help"       => \&usage,
  "version"    => \&version,
  "ref=s"      => \$ref,
  "telomere=s" => \$telomere,
  "debug"      => \$debug,
) or usage(1);
             
$ref or err("Please supply reference genome with --ref");
$ref_index = $ref . ".fai" unless $ref =~ m/\.fai$/;
-r $ref_index or err("Can't see '$ref' index. Run 'samtools faidx $ref' ?"); 
!@ARGV and -t STDIN and err("Please provide or pipe a SAM file to $EXE");

#----------------------------------------------------------------------
# main
msg("$EXE $VERSION by $AUTHOR");

# get a hash of { seqname => length }
msg("Loading: $ref");
my $len = fai_to_dict($ref_index);
msg(Dumper($len)) if $debug;
msg("Found", scalar keys %$len, "sequences in $ref");

my $total=0;
my $addedSequence = {};
my $extendleft=0;
my $extendright=0;
my $repeatsadded=0;
my $addedlen=0;
# at most this number of bases can lack telomeric repeat
# intended for dealing with read errors
my $telosearchspace =  3 * length($telomere); 
my $revcomptelomere=reverse $telomere;
$revcomptelomere =~ tr/ACGTacgt/TGCAtgca/;

# read SAM one line ar a time
while (my $line = <ARGV>) {
  # SAM header
  if ($line =~ m/^@/) {
    # print $line;
    # $header++;
    next;
  }
  $total++;
  my @sam = split m/\t/, $line;
  # ensure there is softclipped alignment before heavyweight parsing
  my $isSoftClipped = ($sam[SAM_CIGAR] =~ /\dS/);
  my $isHardClipped = ($sam[SAM_CIGAR] =~ /\dH/);
  my $isPrimaryAlignment = !($sam[SAM_FLAG] & 2048 or $sam[SAM_FLAG] & 256);
  if ($isSoftClipped and not $isHardClipped and $isPrimaryAlignment ) {
    my $forwardStrand = !($sam[SAM_FLAG] & 16);
    my $start = $sam[SAM_POS];
    my $readlen = length($sam[SAM_SEQ]);
    my $end = $start + $readlen  - 1;
    my $contigname = $sam[SAM_RNAME];
    my $contiglen = $len->{$contigname} or err("Reference", $contigname, "not in '$ref'");
    my ($SL, undef, $SR) 
      = ($sam[SAM_CIGAR] =~ m/ ^ (?:(\d+)S)? (.*?) (?:(\d+)S)? $/x);
    $SL ||= 0; $SR ||= 0;

    ## Add left ends
    if ($forwardStrand and $SL > $start) {
      my $potentialSequence = substr($sam[SAM_SEQ], 0, $SL);
      my $startswithtelomere = ($potentialSequence =~ /^\w{0,\Q$telosearchspace\E}\Q$revcomptelomere/);
      next if not $startswithtelomere;

      # Check whether potentialSequence starts with inverted telomeric repeat?
      if (exists($addedSequence->{$contigname . "L"})) {
        if ($SL > length($addedSequence->{$contigname . "L"})){
          $addedSequence->{$contigname . "L"} = $potentialSequence;
        }
      } else {
        $addedSequence->{$contigname . "L"} = $potentialSequence;
      } 
    } elsif ( ! $forwardStrand ){
    ## Add rightmost ends
      # Check whether alignment ends at the rightmost end of contig
      # adjust end of alignment by indels and deletion length
      # this is obtained by parsing the CIGAR string
      # code adapted from https://github.com/holmeso/adamaperl/blob/56fee71f0c431e3b98ef5f3fc1a2a7213e755e14/lib/QCMG/SamTools/Bam/Alignment.pm#L569
      my $adjust = 0;
      my @cigar  = $sam[SAM_CIGAR] =~ /(\d+)(\w)/g;
      while (@cigar) {
        my ($len,$op) = splice(@cigar,0,2);
        $adjust += $len if $op eq 'I';
        $adjust -= $len if $op eq 'D';
      }
      $adjust += $SL + $SR;
      my $adjustedend = $start + $readlen - $adjust - 1;

      if ( $adjustedend == $contiglen ){
        my $potentialSequence = substr($sam[SAM_SEQ], $readlen - $SR, $readlen + 1);
        my $endswithtelomere = ($potentialSequence =~ /\Q$telomere\E\w{0,\Q$telosearchspace\E}$/);
        next if not $endswithtelomere;
        if (exists($addedSequence->{$contigname . "R"})) {
          if ($SR > length($addedSequence->{$contigname . "R"})){
            $addedSequence->{$contigname . "R"} = $potentialSequence;
          }
        } else {
          $addedSequence->{$contigname . "R"} = $potentialSequence;
        } 
      }
    }
  }
}


my @aux = undef;
my ($name, $seq, $qual);
open my $fh, '<', $ref or err("Can't read FASTA $ref");
while (($name, $seq, $qual) = readfq( $fh, \@aux)) {
  if (exists($addedSequence->{$name . "L"})){
    my $extendingseq = $addedSequence->{$name . "L"};
    $seq = $extendingseq . $seq;
    $repeatsadded += () = $extendingseq =~ /\Q$revcomptelomere/gi;
    $addedlen += length($extendingseq);
    $extendleft++;
  }
  if (exists($addedSequence->{$name . "R"})){
    my $extendingseq = $addedSequence->{$name . "R"};
    $seq .= $extendingseq;
    $repeatsadded += () = $extendingseq =~ /\Q$telomere/gi;
    $addedlen += length($extendingseq);
    $extendright++;
  }
  print ">" . $name . "\n" . $seq . "\n";
}

# stats
msg("Total SAM records $total");
msg("Contigs extended on left side: $extendleft");
msg("Contigs extended on right side: $extendright");
msg("Telomeric repeats added: $repeatsadded");
msg("Total bases added: $addedlen");


msg("Done.");

#----------------------------------------------------------------------
sub fai_to_dict {
  my($fname) = @_;
  open my $fai, '<', $fname or err("Can't read FAI '$fname'");
  my $len = {};
  while (<$fai>) {
    my($name, $bp) = split m/\t/;
    $len->{$name} = $bp;
  }
  close $fai;
  return $len;
}


#----------------------------------------------------------------------
sub version {
  print "$EXE $VERSION\n";
  exit(0);
}

#----------------------------------------------------------------------
sub msg {
  print STDERR "[$EXE] @_\n";
}

#----------------------------------------------------------------------
sub err {
  msg("ERROR:", @_);
  exit(1);
}

#----------------------------------------------------------------------
# Heng Li's readfq: https://github.com/lh3/readfq/blob/master/readfq.pl
sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if !(@$aux);
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}