#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=DRQW1A_2nd_passAligned.out.coverage.bed; data_file_2=temp_long_intron_retention_junction_coordinate.txt; data_file_3=DRQW1A_coverage_temp_long_intron_overlap.txt 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:count_intron_read_long_intron_retention_step1.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %junction; #$junction{junction_id}=[chr][start][end];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") or die "can't open input3 file: $!";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $junction{$array[4]}[0]=$array[0];
      $junction{$array[4]}[1]=$array[1];
      $junction{$array[4]}[2]=$array[2];
}

close IN2;

while(<IN1>) {
      chomp;
      my $ori=$_;
      my @temp=split(/\s+/,$_);
      foreach my $i (keys %junction) {
              if($temp[0] eq $junction{$i}[0]) {
                 if ((($temp[1] <= $junction{$i}[1]) && ($temp[2] >= $junction{$i}[1]) && ($temp[2] <= $junction{$i}[2])) || (($temp[1] >= $junction{$i}[1]) && ($temp[1] <= $junction{$i}[2]) && ($temp[2] >= $junction{$i}[2])) || (($temp[1] >= $junction{$i}[1]) && ($temp[2] <= $junction{$i}[2])) || (($temp[1] <= $junction{$i}[1]) && ($temp[2] >= $junction{$i}[2]))) { 
                   print OUT $ori; print OUT "\n";
                   last;
                 }
              }
     }
}
close IN1;   
