#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use File::Basename;

my ($sam,$eccDNA,$fasta);
my $mapq = 20;
my $overhang = 10;
my $dist = 500;
GetOptions(
  's|sam=s'        => \$sam,
  'e|eccDNA=s'     => \$eccDNA,
  'q|quality=i'    => \$mapq,
  'f|fasta=s'      => \$fasta,
  'o|overhang:i'   => \$overhang,
  'd|dist:i'       => \$dist,
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
}
# chr1    1504182 1504474 7       3       2       SRR6315424.103085058,SRR6315424.52776303

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

my $tracker = 0;
while(<SAM>){
  chomp;
  my @t = split;
  my ($NAME1,$FLAG1,$CHR1,$START1,$MAPQ1,$CIGAR1,$SEQ1) = @t[0,1,2,3,4,5,9];
  my $read2 = <SAM>;
  @t = split("\t",$read2);
  my ($NAME2,$FLAG2,$CHR2,$START2,$MAPQ2,$CIGAR2,$SEQ2) = @t[0,1,2,3,4,5,9];
  $tracker ++;
  print STDERR "$tracker lines have been processed...\n" unless ($tracker%1000000);
# quality filter  
  next if $CHR1 ne $CHR2;
  next if ($MAPQ1 < $mapq || $MAPQ2 < $mapq);
  next if ($CIGAR1 !~ /[SH]/ && $CIGAR2 !~ /[SH]/);
# determine which one is splitted, and which one is contineous
  my ($clip1,$clip2) = (0,0);
  $clip1 = $1 if ($CIGAR1 =~ /(\d+)[SH]/);
  $clip2 = $1 if ($CIGAR2 =~ /(\d+)[SH]/);
  
  if($clip2 > $clip1){
    ($NAME1,$FLAG1,$CHR1,$START1,$MAPQ1,$CIGAR1,$SEQ1) = ($NAME2,$FLAG2,$CHR2,$START2,$MAPQ2,$CIGAR2,$SEQ2);
  }

  my ($len,$mat,$id,$query,$start,$end,$target); 
  if($CIGAR1 =~ /^(\d+)[SH]/){
    $len = $1;
    next if $len < $overhang;
    $mat = 0;
    foreach (keys %{$eccDNA{$CHR1}}){
      $id = $_;
      if($eccDNA{$CHR1}->{$id}->{"leftmost"} < $START1 && $eccDNA{$CHR1}->{$id}->{"leftleast"} + 1 >= $START1){
        $mat = 1; last;
      }
    }
    next if $mat == 0;
    
    my $eccDNAleftmost = join(":",$CHR1,$eccDNA{$CHR1}->{$id}->{"leftmost"});
    foreach (sort {${$left{$eccDNAleftmost}->{'right'}}[$a] <=> ${$left{$eccDNAleftmost}->{'right'}}[$b]} 0..$#{$left{$eccDNAleftmost}->{'right'}}){
      $id = ${$left{$eccDNAleftmost}->{'id'}}[$_];
# length of clipped part > variable regions + $overhang
      next if $len < $overhang + $eccDNA{$CHR1}->{$id}->{"leftleast"} - $eccDNA{$CHR1}->{$id}->{"leftmost"};
      $query = substr($SEQ1,0,$len);
      $end = $eccDNA{$CHR1}->{$id}->{"rightmost"} - ($eccDNA{$CHR1}->{$id}->{"leftleast"} - $START1 + 1);
      $start = $end - $len + 1;
      $target = uc($db->seq($CHR1, $start, $end));
      if ($query eq $target){
        if ($eccDNA{$CHR1}->{$id}->{"leftmost"} < $START2 && $eccDNA{$CHR1}->{$id}->{"rightmost"} >= $START2 + &cigar2glen($CIGAR2) - 1){
          my $STRAND1 = $FLAG1 & 0x10 ? "-" : "+";
          my $STRAND2 = $FLAG2 & 0x10 ? "-" : "+";
          next if $STRAND1 eq $STRAND2;
          if($START1+&cigar2glen($CIGAR1)-1-$eccDNA{$CHR1}->{$id}->{"leftmost"} + $eccDNA{$CHR1}->{$id}->{"rightleast"}-$START2+1 <= $dist){
            push @{$OUT{$id}},$NAME1;
            last;
          }
        }
      }
    }
  }elsif($CIGAR1 =~ /(\d+)[HS]$/){
    $len = $1;
    next if $len < $overhang;
    $mat = 0;
    foreach(keys %{$eccDNA{$CHR1}}){
      $id = $_;
      if($eccDNA{$CHR1}->{$id}->{"rightleast"} <= $START1 + &cigar2glen($CIGAR1) - 1  && $eccDNA{$CHR1}->{$id}->{"rightmost"} >= $START1 + &cigar2glen($CIGAR1) - 1){
         $mat = 1; last;
      }
    }
    next if $mat == 0;

    my $eccDNArightleast = join(":",$CHR1,$eccDNA{$CHR1}->{$id}->{"rightleast"});
    foreach (sort {${$right{$eccDNArightleast}->{'left'}}[$b] <=> ${$right{$eccDNArightleast}->{'left'}}[$a]} 0..$#{$right{$eccDNArightleast}->{'left'}}){
      $id = ${$right{$eccDNArightleast}->{'id'}}[$_];
      next if $len < $overhang + $eccDNA{$CHR1}->{$id}->{"leftleast"} - $eccDNA{$CHR1}->{$id}->{"leftmost"};
      $query = substr($SEQ1,-$len);
      $start = $START1 + &cigar2glen($CIGAR1) - 1 - $eccDNA{$CHR1}->{$id}->{"rightleast"} + $eccDNA{$CHR1}->{$id}->{"leftmost"} + 1;
      $end = $start + $len - 1;
      $target = uc($db->seq($CHR1, $start, $end));
      if ($query eq $target){
        if ($eccDNA{$CHR1}->{$id}->{"leftmost"} < $START2 && $eccDNA{$CHR1}->{$id}->{"rightmost"} >= $START2 + &cigar2glen($CIGAR2) - 1){
          my $STRAND1 = $FLAG1 & 0x10? "-" : "+";
          my $STRAND2 = $FLAG2 & 0x10? "-" : "+";
          next if $STRAND1 eq $STRAND2;
          if($START2+&cigar2glen($CIGAR2)-1-$eccDNA{$CHR1}->{$id}->{"leftmost"} + $eccDNA{$CHR1}->{$id}->{"rightleast"}-$START1+1 <= $dist){
            push @{$OUT{$id}},$NAME1;
            last;
          }
        }
      }
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
  -s --sam         STR  concordant.pe.bam by bwa-mem [STDIN]
                        1. primary supporting reads for -e should be excluded to avoid double counting
                        2. ordered by read name
                        3. supplementary alignment (hard clipped) should be excluded because H cannot be handled in this version
  -e --eccDNA      STR  eccDNA file
  -q --quality     INT  minimal mapping quality [20]
  -f --fasta       STR  fasta file
  -o --overhang    INT  at least n-bp mapped to non-variable breakpoint region [10]
  -d --dist        INT  max distance between two disc reads [500]
  -h --help             Print this help information
HELP
    exit(-1);
}

