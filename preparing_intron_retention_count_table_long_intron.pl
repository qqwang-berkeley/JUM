#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=count_table_index1_long_intron (generated from R); data_file_2=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples.txt; data_file_3=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_long_intron_retention_left_count.txt; data_file_4=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_long_intron_retention_left_count_list.txt; data_file_5=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_long_intron_retention_right_count.txt; data_file_6=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_long_intron_retention_right_count_list.txt; threshold= #; 

#more than threshold reads both adjacent sections

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:making_intron_retention_count_table_long_intron.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4> <data_file_5> <data_file_6> <threshold>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4, $data_file_5, $data_file_6, $threshold)  = @ARGV;

my $o=0;
my %hash_left; #$hash_left{junction_id}=[chr,strand,start,end,count_no_intron,count_with_intron_span_1]
my %hash_right;#$hash_right{junction_id}=[chr,strand,start,end,count_no_intron,count_with_intron_span_2]
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";
open(OUT2, ">$data_file_4") || die "can't open OUT2 file";
open(OUT3, ">$data_file_5") || die "can't open OUT3 file";
open(OUT4, ">$data_file_6") || die "can't open OUT6 file";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash_left{$array[4]}[0]=$array[0];
      $hash_left{$array[4]}[1]=$array[1];
      $hash_left{$array[4]}[2]=$array[2];
      $hash_left{$array[4]}[3]=$array[3];
      $hash_left{$array[4]}[4]=$array[5];
      $hash_left{$array[4]}[5]=0;
      $hash_right{$array[4]}[0]=$array[0];
      $hash_right{$array[4]}[1]=$array[1];
      $hash_right{$array[4]}[2]=$array[2];
      $hash_right{$array[4]}[3]=$array[3];
      $hash_right{$array[4]}[4]=$array[5];
      $hash_right{$array[4]}[5]=0;
}

close IN2;

while(<IN1>) {
      chomp;
      my @temp=split(/\s+/,$_);
      my @name=split(/-/,$temp[0]);
      if(exists $hash_left{$name[0]}) {
           if($name[1] == 1) {
                  $hash_left{$name[0]}[5]=$temp[1];
           }
           if($name[1] == 2) {
                  $hash_right{$name[0]}[5]=$temp[1];
          }
      }
}
close IN1;

foreach my $n (keys %hash_left) {
         if(($hash_left{$n}[5] >= $threshold) && ($hash_right{$n}[5] >= $threshold)) {
         print OUT "intronretentionleft-"; print OUT $n; print OUT ":"; print OUT "001"; print OUT "\t"; 
         print OUT $hash_left{$n}[4]; print OUT "\n";
         print OUT "intronretentionleft-"; print OUT $n; print OUT ":"; print OUT "002"; print OUT "\t";
         print OUT $hash_left{$n}[5]; print OUT "\n";
         print OUT2 "intronretentionleft-"; print OUT2 $n; print OUT2 "\n";

         print OUT3 "intronretentionright-"; print OUT3 $n; print OUT3 ":"; print OUT3 "001"; print OUT3 "\t";
         print OUT3 $hash_right{$n}[4]; print OUT3 "\n";
         print OUT3 "intronretentionright-"; print OUT3 $n; print OUT3 ":"; print OUT3 "002"; print OUT3 "\t";
         print OUT3 $hash_right{$n}[5]; print OUT3 "\n";
         print OUT4 "intronretentionright-"; print OUT4 $n; print OUT4 "\n";
        }
}


