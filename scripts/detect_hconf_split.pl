#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw (min);

my $sam;
my $mapq = 20;
GetOptions(
  'i|input=s'      => \$sam,
  'q|quality=i'    => \$mapq,
  'h|help'         => sub{usage()}
) || usage();

if(!defined $sam){
  $sam="-";
  open SAM,$sam;
}elsif(-B $sam){
  open SAM, "samtools view $sam |" or die "Can't open $sam: $!";
}else{
  open SAM, $sam or die "Can't open $sam: $!";
}

while(<SAM>){
## READ SPLIT READS
  chomp;
  my @t = split;  
  my ($rname1,$flag1,$chr1,$start1,$mapq1,$cigar1) = @t[0..5];
  my $part2 = (<SAM>);
  @t = split("\t",$part2);
  my ($rname2,$flag2,$chr2,$start2,$mapq2,$cigar2) = @t[0..5];
  my $tlen = &cigar2trlen($cigar1); # get the total read length

## PRE-SELECTION
# different read names
# read name is always ended with _[12] after samblaster
  if($rname1 ne $rname2){print STDERR "Error: split fragments are not ordered based on read names.\n"; exit(-1);}
  $rname1=~s/_[12]$//g; $rname2=~s/_[12]$//g;
  my $read1or2 = $flag1 & 0x80? "2" : "1";
# either primary or supplementary alignment is mapped with low quality
  next if ($mapq1<$mapq || $mapq2<$mapq);
# mapped to different chromosomes
# currently, we are not focusing on circular DNA with different chromosome origins
  next if ($chr1 ne $chr2);
# mapped to different strand (potential inversions)
  my $strand1 = $flag1 & 0x10 ? "-" : "+";
  my $strand2 = $flag2 & 0x10 ? "-" : "+";
  next if $strand1 ne $strand2;
# characterize each alignment
# for eccDNAs of small size, the derived reads may cover the breakpoint for >1 times.
  my ($cigartemp,$lc,$rc);
  my ($status1,$rlen1,$glen1);
  if($cigar1 =~ /^(\d+)[HS].*M(\d+)[HS]$/){
    ($cigartemp,$lc,$rc) = ($cigar1,$1,$2);
    $cigartemp =~ s/\d+[HS]//g;
    if ($lc > $rc){
      $tlen = $tlen - $rc;
      $status1 = "left";
    }elsif($lc < $rc){
      $tlen = $tlen - $lc;
      $status1 = "right";
    }else{
      next;
    }
    $rlen1 = &cigar2rlen($cigartemp);
    $glen1 = &cigar2glen($cigartemp);
  }elsif($cigar1 =~ /^\d+[HS](.*\d+M)$/){
    $status1 = "left";
    $rlen1 = &cigar2rlen($1);
    $glen1 = &cigar2glen($1);
  }elsif($cigar1 =~/^(\d+M.*)\d+[HS]$/){
    $status1 = "right";
    $rlen1 = &cigar2rlen($1);
    $glen1 = &cigar2glen($1);
  }

  my ($status2,$rlen2,$glen2);
  if($cigar2 =~ /^(\d+)[HS].*M(\d+)[HS]$/){
    ($cigartemp,$lc,$rc) = ($cigar2,$1,$2);
    $cigartemp =~ s/\d+[HS]//g;
    if ($lc > $rc){
      $tlen = $tlen - $rc;
      $status2 = "left";
    }elsif($lc < $rc){
      $tlen = $tlen - $lc;
      $status2 = "right";
    }else{
      next;
    }
    $rlen2 = &cigar2rlen($cigartemp);
    $glen2 = &cigar2glen($cigartemp);
  }elsif($cigar2 =~ /^\d+[HS](.*\d+M)$/){
    $status2 = "left";
    $rlen2 = &cigar2rlen($1);
    $glen2 = &cigar2glen($1);
  }elsif($cigar2 =~/^(\d+M.*)\d+[HS]$/){
    $status2 = "right";
    $rlen2 = &cigar2rlen($1);
    $glen2 = &cigar2glen($1);
  }
# two parts are likely not complementary
  next if ($status1 eq $status2);
# (discard) $offset < 0: two parts are not connected
# (keep)    $offset = 0; two parts are seamlessly connected
# (keep)    $offset > 0: homologous regions around the breakpoints
  my $offset = $rlen1 + $rlen2 - $tlen;
  next if ($offset < 0);

## WORKING
# format: name read_index strand chromosome left_least(0-base) left_inside(1-base) right_inside(0-base) right_most(1-base) offset_r offset_g_left offset_g_right
  my ($offset_left,$offset_right) = (0,0);
  my ($base_left,$base_right) = (0,0,0,0);
  my ($offset_left_chars,$offset_right_chars);
  my ($offset_left_comm,$offset_right_comm) = (0,0);
  my ($chars1,$chars2);
  if($status1 eq "left" && $status2 eq "right"){
    if(($cigar1 !~ /[DNI]/ && $cigar2 !~ /[DNI]/) || $offset == 0){
      next if ($start1 - 1 + $offset) >= ($start2 + $glen2 - 1);
      print join("\t", $rname1, $read1or2, $strand1, $chr1);
      print "\t";
      print join("\t",$start1 - 1 + $offset, $start1 + $offset - 1 + $glen1 - $offset);
      print "\t";
      print join("\t",$start2 - 1, $start2 + $glen2 - 1);
    }else{
      $chars1 = &cigar2chars($cigar1);
      for(my $i=0; ;$i++){
        $base_left ++ if $chars1 =~ /^[MI]/;
        $offset_left ++ if $chars1 =~ /^[MDN]/;
        $chars1 =~ s/^[MIDN]//;
        last if $base_left == $offset;
      }
      $offset_left += length($1) if $chars1 =~ /^([DN]+)/;
      $chars2 = &cigar2chars($cigar2);
      for(my $i=0; ;$i++){
        $base_right ++ if $chars2 =~ /[MI]$/;
        $offset_right ++ if $chars2 =~ /[MDN]$/;
        $chars2 =~ s/[MIDN]$//;
        last if $base_right == $offset;
      }  
      $offset_right += length($1) if $chars2 =~ /([DN]+)$/;
      $offset_left_chars = substr(&cigar2chars($cigar1),0,$offset_left);
      $offset_left_comm = $1 if($offset_left_chars =~ /(\d+)M$/);
      $offset_right_chars = substr(&cigar2chars($cigar2),-$offset_left);
      $offset_right_comm = $1 if($offset_right_chars =~ /^(\d+)M/);
      $offset = min($offset_left_comm,$offset_right_comm);
      
      next if ($start1 - 1 + $offset_left) >= ($start2 + $glen2 - 1);
      print join("\t", $rname1, $read1or2, $strand1, $chr1);
      print "\t";
      print join("\t",$start1 - 1 + $offset_left, $start1 + $offset_left - 1 + $glen1 - $offset_left);
      print "\t";
      print join("\t",$start2 - 1, $start2 + $glen2 - 1);
    }
  }else{
    if(($cigar1 !~ /[DNI]/ && $cigar2 !~ /[DNI]/) || $offset == 0){
      next if ($start2 - 1 + $offset) >= ($start1 + $glen1 - 1);
      print join("\t", $rname1, $read1or2, $strand1, $chr1);
      print "\t";
      print join("\t",$start2 - 1 + $offset, $start2 + $offset - 1 + $glen2 - $offset);
      print "\t";
      print join("\t",$start1 - 1, $start1 + $glen1 - 1);
    }else{
      $chars2 = &cigar2chars($cigar2);
      for(my $i=0; ;$i++){
        $base_left ++ if $chars2 =~ /^[MI]/;
        $offset_left ++ if $chars2 =~ /^[MDN]/;
        $chars2 =~ s/^[MIDN]//;
        last if $base_left == $offset;
      }  
      $offset_left += length($1) if $chars2 =~ /^([DN]+)/;
      $chars1 = &cigar2chars($cigar1);
      for(my $i=0; ;$i++){
        $base_right ++ if $chars1 =~ /[MI]$/;
        $offset_right ++ if $chars1 =~ /[MDN]$/;
        $chars1 =~ s/[MIDN]$//;
        last if $base_right == $offset;
      }
      $offset_right += length($1) if $chars1 =~ /([DN]+)$/;
      $offset_left_chars = substr(&cigar2chars($cigar2),0,$offset_left);
      $offset_left_comm = $1 if($offset_left_chars =~ /(\d+)M$/);
      $offset_right_chars = substr(&cigar2chars($cigar1),-$offset_left);
      $offset_right_comm = $1 if($offset_right_chars =~ /^(\d+)M/);
      $offset = min($offset_left_comm,$offset_right_comm);

      next if ($start2 - 1 + $offset_left) >= ($start1 + $glen1 - 1);
      print join("\t", $rname1, $read1or2, $strand1, $chr1);
      print "\t";
      print join("\t",$start2 - 1 + $offset_left, $start2 + $offset_left - 1 + $glen2 - $offset_left);
      print "\t";
      print join("\t",$start1 - 1, $start1 + $glen1 - 1);
    }
  }
  print "\t$offset\n";
}

## SUBROUTINES
# infer the length of the original read
sub cigar2trlen{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MISH]/g);
    $length;
}
# length of mappable part of a read
sub cigar2rlen{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MI]/g);
    $length;
}
# length of mappable region
sub cigar2glen{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MDN]/g);
    $length;
}
# cigar to chars
sub cigar2chars {
  my ($cigar) = @_;
  my $char;
  $cigar =~ s/\d+[SH]//g;
  foreach(;$cigar =~ s/^(\d+)([MDIN])//;){
    $char .= $2 x $1;
  }
  $char;
}
# help message
sub usage{
my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -i .split.[bs]am > OUTPUT
  If INPUT isn't specified, input from STDIN
Options:
  -i --input       STR  split reads in .[sb]am format (eg., split reads after bwa-mem | samblaster)
  -q --quality     INT  minimal mapping quality [20]
  -h --help             Print this help information
Release Notes:
  \@version 1.0
    beta version
HELP
    exit(-1);
}

