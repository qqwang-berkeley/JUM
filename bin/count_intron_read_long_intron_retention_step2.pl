#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=DRQW1A_coverage_temp_long_intron_overlap.txt; data_file_2=temp_long_intron_retention_junction_coordinate.txt; data_file_3=temp_long_intron_retention_junction_coordinate_with_read_num.txt 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:count_intron_read_long_intron_retention_step2.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %section; #$section{chr}{start}{end}=[#num];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") or die "can't open input3 file: $!";

my $num;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_); 
         if($array[2]-$array[1] > 1) {
              for($num=0; $num<=$array[2]-$array[1]-1; $num++) {
                   $section{$array[0]}{$array[1]+$num}{$array[1]+$num+1}[0]=$array[3];
              }
         }
         else {
                   $section{$array[0]}{$array[1]}{$array[2]}[0]=$array[3];
         }
}

close IN1;

while(<IN2>) {
      chomp;
      my @temp=split(/\s+/,$_); 
      if($temp[2]-$temp[1] > 0) {
           if($temp[2]-$temp[1] > 1) {
              for($num=0; $num<=$temp[2]-$temp[1]-1; $num++) {
                   #if(exists $section{$temp[0]}{$temp[1]+$num}{$temp[1]+$num+1}) {
                       print OUT $temp[0]; print OUT "\t"; print OUT $temp[1]+$num; print OUT "\t"; print OUT $temp[1]+$num+1; print OUT "\t"; print OUT $temp[3]; print OUT "\t"; print OUT $temp[4]; print OUT "\t"; print OUT $section{$temp[0]}{$temp[1]+$num}{$temp[1]+$num+1}[0]; print OUT "\n";
              }
         }
          else {
                       print OUT $temp[0]; print OUT "\t"; print OUT $temp[1]; print OUT "\t"; print OUT $temp[2]; print OUT "\t"; print OUT $temp[3]; print OUT "\t"; print OUT $temp[4]; print OUT "\t"; print OUT $section{$temp[0]}{$temp[1]}{$temp[2]}[0]; print OUT "\n";
        }
     }
}

close IN2; 
