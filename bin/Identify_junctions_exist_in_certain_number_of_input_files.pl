#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=UNION_junc_coor_with_junction_ID_more_than_5_reads.txt; data_file_2=DRQW1A_2nd_pass_junction_counts.txt; data_file_3=DRQW1B_2nd_pass_junction_counts.txt; data_file_4=DRQW1C_2nd_pass_junction_counts.txt ; data_file_4=...; data_file_5=...; data_file_n-1=OUTPUT; data_file_n-1=$file_number; data_file_n=$threshold;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:Identify_junctions_exist_in_certain_number_of_input_files.pl <input_1> <input_2> ... <output_file> <file_number> <threshold>\n";
}

my $mark;
my $loop;
my $o=1;
my %hash; #$hash{junction_ID}=[$chr,$strand, $start, $end, $index];
my @file;
$file[0]=0;
my $filehandle;

for($mark=1;$mark<=$len-3;$mark++) {
        open($file[$mark], "<$ARGV[$mark-1]") || die "can't open intput file $mark";
}

open(OUT, ">$ARGV[$len-3]") || die "can't open OUT file";

my $file_num = $ARGV[$len-2];
my $threshold = $ARGV[$len-1];

$filehandle=$file[1];

while(<$filehandle>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[4]}[0]=0;
}

close $filehandle;

       for ($loop=2; $loop<=$len-3; $loop++) {
             $filehandle = $file[$loop];
             while(<$filehandle>) {
                     chomp;
                     my @Xtemp=split(/\s+/,$_);
                     if((exists $hash{$Xtemp[4]}) && ($Xtemp[5] >= $threshold))  {
                             $hash{$Xtemp[4]}[0]=$hash{$Xtemp[4]}[0]+1;
                     }
            }
             close $filehandle;
      }


print $len-4; print " "; print "input samples to consider"; print "\n";
foreach my $n (keys %hash) {
             if ($hash{$n}[0] >= $file_num) {
                          print OUT $n; print OUT "\n";
            }
}
