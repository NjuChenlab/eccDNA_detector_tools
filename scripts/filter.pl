#!/usr/bin/perl
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($originalFile);
my ($originalFields, $targetFields, $mode, $sep)=("1", "1", "e", "");

sub usage{
  my $scriptName = basename $0;
print <<HELP;
Usage:
  perl $scriptName -o originFile.tsv -1 1,4 -m i|include targetFile.tsv >filtered.tsv
  if targetFile.tsv isn't specified, input is from STDIN, output to STDOUT
Options:
  -o --originFile     The original file containing fields (specified by --originFields) used to include or exclude lines in targetFile.tab
  -1 --originFields   Comma-separated field list specifying which fileds in the originFile.tab to be used to include or exclude, 1-based start [1]
                      The element of the list can be a single column number or a range with nonnumeric char as separator
                      To specify the last column left the range right margin blank
                      If continuous range specified like '1-3-6', the first range '1-3' will be output
                      eg.:
                        -1 1,4          output columns 1,4
                        -1 1-4,6..8     output columns 1,2,3,4,6,7,8
                        -1 1,4,6-       output columns 1,4,6,7,... last column
                        -1 1-3-6        output columns 1,2,3
  -2 --targetFields   Comma-separated field list specifying which fileds in the targetFile.tab are used to include or exclude lines, 1-based start [1]
                      More description about --targetFields, see --originFields
  -m --mode           To include or exclude lines in targetFile.tab, it can be i|include or e|exclude[e]
  -s --separator      (Optional)A separator to join the fields specified, if necessary[Empty string]
  -h --help           Print this help information screen
HELP
    exit(-1);
}

GetOptions(
  'o|originFile=s'    => \$originalFile,
  '1|originFields=s'  => \$originalFields,
  '2|targetFields=s'  => \$targetFields,
  'm|mode=s'          => \$mode,
  's|separator=s'     => \$sep,
  'h|help'            => sub{&usage()}
)||usage();

die "Plese specify -o to offer a file containing the list your want to include or exclude!\n" if !defined $originalFile;
$ARGV[0]='-' unless defined $ARGV[0];
open ORI,"$originalFile" or die "$originalFile:$!";
open TAR,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

my %originalTag;
while(<ORI>){
  chomp;
  my @fields = split "\t";
  my $tag = join $sep, @{&arraySplice(\@fields, $originalFields)->{array}};
  $originalTag{$tag} = '';
}

if($mode =~ /i|include/){
  while(<TAR>){
    chomp;
    my @fields = split "\t";
    my $tag = join $sep, @{&arraySplice(\@fields,$targetFields)->{array}};
    say if( defined $originalTag{$tag} );
  }
}elsif($mode =~ /e|exclude/){
  while(<TAR>){
    chomp;
    my @fields = split "\t";
    my $tag = join $sep, @{&arraySplice(\@fields,$targetFields)->{array}};
    say if( !defined $originalTag{$tag} );
  }
}else{
  die "Please specify the correct mode(i|include|e|exclude)\n";
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

