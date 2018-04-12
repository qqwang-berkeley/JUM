#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=UNION_junc_coor_with_junction_ID_more_than_5_reads.txt; data_file_2=DRQW1A_2nd_pass_junction_counts.txt; data_file_3=DRQW1B_2nd_pass_junction_counts.txt; data_file_4=DRQW1C_2nd_pass_junction_counts.txt ; data_file_4=...; data_file_5=...; data_file_6=OUTPUT 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:Identify_junctions_exist_in_all_samples.pl <input_1> <input_2> ... <output_file>\n";
}

my $mark;
my $loop;
my $o=1;
my %hash; #$hash{junction_ID}=[$chr,$strand, $start, $end, $index];
my @file;
$file[0]=0;
my $filehandle;

for($mark=1;$mark<=$len-1;$mark++) {
        open($file[$mark], "<$ARGV[$mark-1]") || die "can't open intput file $mark";
}

open(OUT, ">$ARGV[$len-1]") || die "can't open OUT file";

$filehandle=$file[1];

while(<$filehandle>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[4]}[0]=$array[0];
      $hash{$array[4]}[1]=$array[3];
      $hash{$array[4]}[2]=$array[1];
      $hash{$array[4]}[3]=$array[2];
      $hash{$array[4]}[4]=0;
}

close $filehandle;

       for ($loop=2; $loop<=$len-1; $loop++) {
             $filehandle = $file[$loop];
             while(<$filehandle>) {
                     chomp;
                     my @Xtemp=split(/\s+/,$_);
                     if(exists $hash{$Xtemp[4]}) {
                             $hash{$Xtemp[4]}[4]=$hash{$Xtemp[4]}[4]+1;
                     }
            }
             close $filehandle;
      }


print $len-2; print " "; print "input samples for JUM analyses"; print "\n";
foreach my $n (keys %hash) {
             if ($hash{$n}[4] == $len-2) {
                          print OUT $n; print OUT "\n";
            }
}
