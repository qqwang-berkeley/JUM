#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );

#data_file_1=left_span_long_intron_retention_AS_differential_combined_pvalue_0_05.txt; data_file_2=right_span_long_intron_retention_AS_differential_combined_pvalue_0_05.txt; data_file_3=significant_long_intron_retention_splicing_event.txt; data_file_4=temp_long_intron_list.txt

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:prepare_significant_long_intron_retention_event_after_DEXSeq_analyses.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4>\n";
}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4)  = @ARGV;

my %hash_left; #$hash{junction_id}{sub_junction_id}=[left_whole_line,left_delta_change,right_whole_line,right_delta_change,indicator_for_right_matching_left,indicator_for_change_equivalent_in_both_left_and_right];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file: $!";
open(OUT1, ">$data_file_4") || die "can't open OUT1 file: $!";

while(<IN1>) {
      chomp;
      my $origin=$_;
      my @left=split(/\s+/,$origin);
      my @leftside=split(/-/,$left[0]);
      my @leftid=split(/:/,$leftside[1]);
      if($leftid[1] eq "E002") {
        $hash_left{$leftid[0]}[0]=$origin;
        $hash_left{$leftid[0]}[1]=$left[10];
        $hash_left{$leftid[0]}[2]=0;
        $hash_left{$leftid[0]}[3]=0;
      }
}
close IN1;

while(<IN2>) {
      chomp;
      my $rightori=$_;
      my @right=split(/\s+/,$rightori);
      my @rightside=split(/-/,$right[0]);
      my @rightid=split(/:/,$rightside[1]);
      if($rightid[1] eq "E002") {
          #if((exists $hash_left{$rightid[0]}) && (abs($right[10]-$hash_left{$rightid[0]}[1])/max(abs($right[10]),abs($hash_left{$rightid[0]}[1])) <= 1)) {
           if(exists $hash_left{$rightid[0]}) {
              $hash_left{$rightid[0]}[2]=$rightori;
              $hash_left{$rightid[0]}[3]=1;
          }
     }
}

close IN2;

foreach my $ke (keys %hash_left) {
                   if($hash_left{$ke}[3]==1) {
                   print OUT $ke; print OUT "\n";
                   print OUT "subevent"; print OUT "\t"; print OUT $hash_left{$ke}[0]; print OUT "\n"; print OUT "subevent"; print OUT "\t"; print OUT $hash_left{$ke}[2]; print OUT "\n";
                   print OUT1 $ke; print OUT1 "\n";
            }
}
                 

