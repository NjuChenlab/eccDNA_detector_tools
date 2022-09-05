#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use File::Basename;

my ($sam,$eccDNA,$fasta);
my $mapq = 20;
my $overhang = 10;
GetOptions(
  's|sam=s'        => \$sam,
  'e|eccDNA=s'     => \$eccDNA,
  'q|quality=i'    => \$mapq,
  'f|fasta=s'      => \$fasta,
  'o|overhang:i'   => \$overhang,
  'h|help'         => sub{usage()}
) || usage();

if(!defined $sam){
  $sam="-";
  open SAM,$sam;
}elsif(-B $sam){
  open SAM, "samtools view $sam|" or die "Can't open $sam: $!";
}else{
  open SAM, "$sam" or die "Can't open $sam: $!";
}

if(!defined $eccDNA){
  print STDERR "Pls give eccDNA files.\n";
  exit(-1);
}else{
  open DNA, $eccDNA;
# chr1    1504182 1504474 7       3       2       SRR6315424.103085058,SRR6315424.52776303
}

my $db;
if(!defined $fasta){
  print STDERR "Pls give .fa format reference genome.\n";
  exit(-1);
}else{
  $db = Bio::DB::Fasta->new($fasta);
}

my %eccDNA;
my %left;
my %right;
my %OUT;
while(<DNA>){
  chomp;
  my @t = split;
  my ($chr,$leftmost,$rightleast,$name,$offset) = @t[0,1,2,3,4];
  my $leftleast = $leftmost + $offset;
  my $rightmost = $rightleast + $offset; 
   
  $eccDNA{$chr}->{$name}->{"leftmost"}   = $leftmost; 
  $eccDNA{$chr}->{$name}->{"leftleast"}  = $leftleast; 
  $eccDNA{$chr}->{$name}->{"rightleast"} = $rightleast; 
  $eccDNA{$chr}->{$name}->{"rightmost"}  = $rightmost; 
  
  my $leftid = join(":",$chr,$leftmost);
  push @{$left{$leftid}->{'right'}},$rightleast;
  push @{$left{$leftid}->{'id'}},$name;

  my $rightid = join(":",$chr,$rightleast);
  push @{$right{$rightid}->{'left'}},$leftmost;
  push @{$right{$rightid}->{'id'}},$name;
  
  @{$OUT{$name}} = ();
}

my ($NAME,$FLAG,$CHR,$START,$MAPQ,$CIGAR,$SEQ);
my $tracker = 0;
while(<SAM>){
  chomp;
  my @t = split;
  ($NAME,$FLAG,$CHR,$START,$MAPQ,$CIGAR,$SEQ) = @t[0,1,2,3,4,5,9];
  $tracker ++;
  print STDERR "$tracker lines have been processed...\n" unless ($tracker%1000000);
# quality filter
  next if $MAPQ < $mapq;
# length of clipped part
  next if $CIGAR !~ /[SH]/;
  my ($len,$mat,$id,$query,$start,$end,$target); 
  if($CIGAR =~ /^(\d+)[SH]/){
    $len = $1;
    next if $len < $overhang;
    $mat = 0;
    foreach (keys %{$eccDNA{$CHR}}){
      $id = $_;
      if($eccDNA{$CHR}->{$id}->{"leftmost"} < $START && $eccDNA{$CHR}->{$id}->{"leftleast"} + 1 >= $START){
        $mat = 1; last;
      }
    }
    next if $mat == 0;

    my $eccDNAleftmost = join(":",$CHR,$eccDNA{$CHR}->{$id}->{"leftmost"});
    foreach (sort {${$left{$eccDNAleftmost}->{'right'}}[$a] <=> ${$left{$eccDNAleftmost}->{'right'}}[$b]} 0..$#{$left{$eccDNAleftmost}->{'right'}}){
  # foreach (0..$#{$left{$eccDNAleftmost}->{'right'}}){
      $id = ${$left{$eccDNAleftmost}->{'id'}}[$_];
# length of clipped part > variable regions + $overhang
      next if $len < $overhang + $eccDNA{$CHR}->{$id}->{"leftleast"} - $eccDNA{$CHR}->{$id}->{"leftmost"};
      $query = substr($SEQ,0,$len);
      $end = $eccDNA{$CHR}->{$id}->{"rightmost"} - ($eccDNA{$CHR}->{$id}->{"leftleast"} - $START + 1);
      $start = $end - $len + 1;
      $target = uc($db->seq($CHR, $start, $end));
      push @{$OUT{$id}},$NAME if $query eq $target;
      last;
    }
  }elsif($CIGAR =~ /(\d+)[HS]$/){
    $len = $1;
    next if $len < $overhang;
    $mat = 0;
    foreach(keys %{$eccDNA{$CHR}}){
      $id = $_;
      if($eccDNA{$CHR}->{$id}->{"rightleast"} <= $START + &cigar2glen($CIGAR) - 1 && $eccDNA{$CHR}->{$id}->{"rightmost"} >= $START + &cigar2glen($CIGAR) - 1){
         $mat = 1; last;
      }
    }
    next if $mat == 0;

    my $eccDNArightleast = join(":",$CHR,$eccDNA{$CHR}->{$id}->{"rightleast"});
    foreach (sort {${$right{$eccDNArightleast}->{'left'}}[$b] <=> ${$right{$eccDNArightleast}->{'left'}}[$a]} 0..$#{$right{$eccDNArightleast}->{'left'}}){
  # foreach (0..$#{$right{$eccDNArightleast}->{'left'}}){
      $id = ${$right{$eccDNArightleast}->{'id'}}[$_];
      next if $len < $overhang + $eccDNA{$CHR}->{$id}->{"leftleast"} - $eccDNA{$CHR}->{$id}->{"leftmost"};
      $query = substr($SEQ,-$len);
      $start = $START + &cigar2glen($CIGAR) - 1 - $eccDNA{$CHR}->{$id}->{"rightleast"} + $eccDNA{$CHR}->{$id}->{"leftmost"} + 1;
      $end = $start + $len - 1;
      $target = uc($db->seq($CHR, $start, $end));
      push @{$OUT{$id}},$NAME if $query eq $target;
      last;
    }
  }
}

foreach(keys %OUT){
  print "$_\t";
  if(scalar @{$OUT{$_}} == 0){
    print "0\t-";
  }else{
    print scalar @{$OUT{$_}};
    print "\t";
    print join(",",@{$OUT{$_}});
  }
  print "\n";
}

## SUBROUTINES
sub cigar2glen{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MDN]/g);
    $length;
}
sub usage{
my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -i .split.[bs]am > OUTPUT
Options:
  -s --sam         STR  singleton.bam by bwa-mem [STDIN]
                        1. primary supporting reads for -e should be excluded to avoid double counting
                        2. ordered by read name
                        3. supplementary alignment (hard clipped) should be excluded
  -e --eccDNA      STR  eccDNA file
  -q --quality     INT  minimal mapping quality [20]
  -f --fasta       STR  fasta file
  -o --overhang    INT  at least n-bp mapped to non-variable breakpoint region [10]
  -h --help             Print this help information

\@release note
  last modified: 05/15/2020
  \@ver 1.1
    some eccDNAs may share one break end, each of these eccDNAs are checked.
    if ties occur and you choose to believe the shortest eccDNA, you should sort -k 1,1 -k 2,2nr -k 3,3n .eccDNA file;
    if ties occur and you choose to believe the longest  eccDNA, you should sort -k 1,1 -k 2,2n -k 3,3nr .eccDNA file;
  \@ver 1.0
    beta version

HELP
    exit(-1);
}

