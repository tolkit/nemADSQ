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
my $minlen    = 0;
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
    "  % minimap2 -a ref.fa telomeric_reads.fa | extract_softlipped_non_tel_seqs.pl --minlen 1000 > soft_clipped_seqs.fasta\n",
    "OPTIONS\n",
    "  --help           This help\n",
    "  --version        Print version and exit\n",
    "  --minlen FASTA      Reference genome - needs FASTA.fai index\n",
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
  "minlen=i"      => \$minlen,
  "telomere=s" => \$telomere,
  "debug"      => \$debug,
) or usage(1);

!@ARGV and -t STDIN and err("Please provide or pipe a SAM file to $EXE");

#----------------------------------------------------------------------
# main
msg("$EXE $VERSION by $AUTHOR");


# at most this number of bases can lack telomeric repeat
# intended for dealing with read errors
my $total=0;
my $softclip_sequences = {};
my $record_count = 0;
my $telosearchspace =  3 * length($telomere); 
my $revcomptelomere=reverse $telomere;
$revcomptelomere =~ tr/ACGTacgt/TGCAtgca/;


# read SAM one line at a time
while (my $line = <ARGV>) {
  # SAM header
  if ($line =~ m/^@/) {
    # print $line;
    # $header++;
    next;
  }
  $total++;
  my @sam = split m/\t/, $line;
  #print "hasdasdlklkmlklo\n" . "$sam[0]" . "\t" . $sam[SAM_CIGAR] . "\n";
  # ensure there is softclipped alignment before heavyweight parsing
  my $isSoftClipped = ($sam[SAM_CIGAR] =~ /\dS/);
  my $isHardClipped = ($sam[SAM_CIGAR] =~ /\dH/);
  my $isPrimaryAlignment = !($sam[SAM_FLAG] & 2048 or $sam[SAM_FLAG] & 256);
  if ($isSoftClipped and not $isHardClipped and $isPrimaryAlignment ) {
    my $readname = $sam[0];
    my $forwardStrand = !($sam[SAM_FLAG] & 16);
    my $start = $sam[SAM_POS];
    my $readlen = length($sam[SAM_SEQ]);
    my $end = $start + $readlen  - 1;
    my ($SL, undef, $SR) 
      = ($sam[SAM_CIGAR] =~ m/ ^ (?:(\d+)S)? (.*?) (?:(\d+)S)? $/x);
    $SL ||= 0; $SR ||= 0;


    next if ($SR < $minlen and $SL < $minlen);
    # print "hasdasdlklkmlklo\n" . "$SL" . "\n" . "$SR" . "\n";

    ## Add left ends
    if ($forwardStrand and $SL > $SR) {
      # print "hasdasdlklkmlklo\n" . "$SL" . "\n" . "$SR" . "\n";
      my $potentialSequence = substr($sam[SAM_SEQ], 0, $SL + 1);
      my $startswithtelomere = ($potentialSequence =~ /^\w{0,\Q$telosearchspace\E}\Q$revcomptelomere/);
      #print ">" . $readname . "\n" . $potentialSequence . "\n";
      next if not $startswithtelomere;
      print ">" . $readname . "\n" . $potentialSequence . "\n";


    } elsif ( ! $forwardStrand ){
      my $potentialSequence = substr($sam[SAM_SEQ], -$SR);
      my $endswithtelomere = ($potentialSequence =~ /\Q$telomere\E\w{0,\Q$telosearchspace\E}$/);
      #print ">" . $readname . "\n" . $potentialSequence . "\n";
      next if not $endswithtelomere;
      print ">" . $readname . "\n" . $potentialSequence . "\n";
    }
  }
}



msg("Done.");




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