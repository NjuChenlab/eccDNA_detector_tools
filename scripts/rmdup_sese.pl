#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw (max);

my $sam;
GetOptions(
  'i|input=s'      => \$sam,
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

my ($name,$md,$string,$mapq,$strand);
my %rmdup;

while(<SAM>){
  chomp;
  next if /^@/;
  my @t = split;
  next if $t[1] & 0x4;
  if(!defined $name){
    $name = $t[0];
    $md = $1 if /\tMD:Z:([\^0123456789ATCG]*)/;
    $t[1] & 0x10 ? $strand = "-" : $strand = "+";
    $string = join(':',@t[2,3,5],$strand,$md);
    $mapq = $t[4];
    $rmdup{$string}->{'mapq'} = $mapq;
    $rmdup{$string}->{'name'} = $name;
  }elsif($name eq $t[0]){
    $md = $1 if /\tMD:Z:([\^0123456789ATCG]*)/;
    $t[1] & 0x10 ? $strand = "-" : $strand = "+";
    $string = $string . ":" . join(':',@t[2,3,5],$strand,$md);
    $mapq = ($mapq + $t[4]) / 2;
    if(exists $rmdup{$string}){
      if ($rmdup{$string}->{'mapq'} < $mapq){
        $rmdup{$string}->{'mapq'} = $mapq;
        $rmdup{$string}->{'name'} = $name;
      }
    }else{
      $rmdup{$string}->{'mapq'} = $mapq;
      $rmdup{$string}->{'name'} = $name;
    }
  }else{
    $name = $t[0];
    $md = $1 if /\tMD:Z:([\^0123456789ATCG]*)/;
    $t[1] & 0x10 ? $strand = "-" : $strand = "+";
    $string = join(':',@t[2,3,5],$strand,$md);
    $mapq = $t[4];
    if(exists $rmdup{$string}){
      if ($rmdup{$string}->{'mapq'} < $mapq){
        $rmdup{$string}->{'mapq'} = $mapq;
        $rmdup{$string}->{'name'} = $name;
      }
    }else{
      $rmdup{$string}->{'mapq'} = $mapq;
      $rmdup{$string}->{'name'} = $name;
    }
  }
}
close SAM;

if(!defined $sam){
  $sam="-";
  open SAM,$sam;
}elsif(-B $sam){
  open SAM, "samtools view -h $sam |" or die "Can't open $sam: $!";
}else{
  open SAM, $sam or die "Can't open $sam: $!";
}

my %ids;
foreach(keys %rmdup){
  my $id = $rmdup{$_}->{'name'};
  delete $rmdup{$_};
  $ids{$id} = 1;
}
while(<SAM>){
  chomp;
  if(/^@/){print "$_\n"; next;}
  my @t = split;
  if(exists $ids{$t[0]}){
    print join("\t",@t);
    print "\n";
  }
}

sub usage{
my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -i .se.ori/nsrt.bam > OUTPUT
  If INPUT isn't specified, input from STDIN
Options:
  -i --input       STR  split reads in .[sb]am format (eg., split reads after bwa-mem | samblaster)
  -h --help             Print this help information
Release Notes:
  \@version 1.0
    rm duplicate for se.s/bam; secondary alignments are taken into consideration 
HELP
    exit(-1);
}
