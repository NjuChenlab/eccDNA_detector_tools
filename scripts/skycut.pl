#!/usr/bin/perl
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
my $fields = "";

GetOptions(
  'f|fields=s'  => \$fields,
  'h|help'      => sub{usage()}
) || usage();

$ARGV[0]='-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

die "please specify the correct -f or --fields\n" if $fields !~ /(\d+\D*\d*,?)+/;

while(<IN>){
  my @fields = split "\t";
  chomp $fields[-1];
  say join "\t", @{&arraySplice(\@fields, $fields)->{array}};
}

sub arraySplice{
  my ($arrayI,$splices)=@_;
  my @splicesA=split ",",$splices;
  my @indexs;
  for my $splice (@splicesA){
    if($splice=~/\D+/){#-f 1- or -7 or 1-4 or -f 4-1 or -f1-1
      if($splice=~/\D+$/){#-f 1-
        my $from=(split /\D+/,$splice)[0];
        die "'$splice' in the fields your specify ($splices) isn't in correct form" unless defined $from;
        push @indexs,$_ for( $from-1..$#$arrayI );
      }elsif( $splice=~/^\D+/ ){#-f -7
        die "'$splice' in the fields your specify ($splices) isn't in correct form" if split /\D+/,$splice !=2;
        my $from=(split /\D+/,$splice)[1];
        push @indexs,$_ for reverse( $from-1..$#$arrayI );
      }else{#-f 1-4 or -f 4-1 or -f 1-1
        my ($from,$to)=split /\D+/,$splice;
        if($from<$to){#-f 1-4
          if($to>@$arrayI){
            say STDERR "Warnning: no $to columns in: ".join "\t",@$arrayI;
            $to=@$arrayI;
          }
          push @indexs,$_ for( ($from-1)..($to-1) );
        }else{#-f 4-1 or -f 1-1
          ($from,$to)=($to,$from);
          if($to>@$arrayI){
            say STDERR "Warnning: no $to columns in: ".join "\t",@$arrayI;
            $to=$#$arrayI;
          }
          push @indexs,$_ for reverse( ($from-1)..($to-1) );
        }
      }
    }else{#-f 1
      if($splice>@$arrayI){
        say STDERR "Warnning: no $splice columns in: ".join "\t",@$arrayI;
      }else{
        push @indexs,($splice-1);
      }
    }
  }
  my @arraySpliced=map{ $arrayI->[$_] }@indexs;
  return {
    "index" => \@indexs,
    "array" => \@arraySpliced
  };
}

sub usage{
  my $scriptName = basename $0;
print <<HELP;
Usage:perl $scriptName INPUT.tsv >OUTPUT.tsv
  if INPUT not specified,input from STDIN, output from STDOUT
  -f --fields     Comma-separated list specifying the fields to output
                  The element of the list can be a single column number or a range with nonnumeric char as separator
                  To specify the last column left the range right margin blank
                  If continuous range specified like '1-3-6', the first range '1-3' will be output
                  eg.:
                      -f 1,4          output columns 1,4
                      -f 1-4,6..8     output columns 1,2,3,4,6,7,8
                      -f 1,4,6-       output columns 1,4,6,7,... last column
                      -f 1-3-6        output columns 1,2,3
  -h --help       Print this help information
HELP
  exit(-1);
}



