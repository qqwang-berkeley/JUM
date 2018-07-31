#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=DRQW1A_2nd_passAligned.out.coverage_with_ID.bed; data_file_2=temp_long_intron_retention_sliding_windowns.txt; data_file_3=temp_long_intron_retention_sliding_windowns_with_read.txt 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:count_intron_read_long_intron_retention.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %section; #$section{section_id}=[chr][start][end][#num];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(IN3, ">$data_file_3") or die "can't open input3 file: $!";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      $section{$array[4]}[0]=$array[0];
      $section{$array[4]}[1]=$array[1];
      $section{$array[4]}[2]=$array[2];
      $section{$array[4]}[3]=$array[3];
}

close IN1;

while(<IN2>) {
      chomp;
      my $ori=$_;
      my @temp=split(/\s+/,$_);
      my $mark="*";
      foreach my $i (keys %section) {
          if(($temp[0] eq $section{$i}[0]) && ($temp[1] >= $section{$i}[1]) && ($temp[2] <= $section{$i}[2])) {
                  $mark=$section{$i}[3];
                  last;
          }
      }
     print IN3 $ori; print IN3 "\t"; print IN3 $mark; print IN3 "\n";
}
close IN2;   
