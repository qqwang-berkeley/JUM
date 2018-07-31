#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=Index1_Aligned.out.spanning_junction_reads_bed12.bed; data_file_2=output; readlength= #

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profile_overhang_from_STAR_output_step_1_no_strand_considered.pl <data_file_1> <data_file_2> <readlength>\n";

}

my ($data_file_1, $data_file_2, $readlength)  = @ARGV;

my %hash; #$hash{chr}{start}{end}=[leftoverhang,rightoverhang];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") || die "can't open OUT file";
my $i;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      my @subarray_lag=split(/,/,$array[10]);
      my @subarray_int=split(/,/,$array[11]);
      my $end=$array[1];
      for($i=0;$i<=$#subarray_lag-1;$i++) {
          my $start_jun=$end+$subarray_lag[$i]+1;
          my $end_jun=$array[1]+$subarray_int[$i+1];
              if (exists $hash{$array[0]}{$start_jun}{$end_jun}) {
                  $hash{$array[0]}{$start_jun}{$end_jun}[0]=$hash{$array[0]}{$start_jun}{$end_jun}[0]+1;
                   if(($subarray_lag[$i] > $hash{$array[0]}{$start_jun}{$end_jun}[1]) && ($subarray_lag[$i] < $readlength)) {
                     $hash{$array[0]}{$start_jun}{$end_jun}[1]=$subarray_lag[$i];
                  }
                   if(($subarray_lag[$i+1] > $hash{$array[0]}{$start_jun}{$end_jun}[2]) && ($subarray_lag[$i+1] < $readlength)) {
                     $hash{$array[0]}{$start_jun}{$end_jun}[2]=$subarray_lag[$i+1];
                  }
              }
              else {
                 $hash{$array[0]}{$start_jun}{$end_jun}[0]=1;
                 if($subarray_lag[$i] < $readlength) {
                 $hash{$array[0]}{$start_jun}{$end_jun}[1]=$subarray_lag[$i];
                 }
                 else {
                 $hash{$array[0]}{$start_jun}{$end_jun}[1]=$readlength-1;
                 }
                 if($subarray_lag[$i+1] < $readlength) {
                 $hash{$array[0]}{$start_jun}{$end_jun}[2]=$subarray_lag[$i+1];
                 }
                 else {
                 $hash{$array[0]}{$start_jun}{$end_jun}[2]=$readlength-1;
                 }
              }
          $end=$end_jun;
      }
}


foreach my $chr (keys %hash) {
            foreach my $start (keys %{$hash{$chr}}) {
                 foreach my $end (keys %{$hash{$chr}{$start}}) {
                      print OUT $chr; print OUT "\t"; print OUT $start; print OUT "\t"; print OUT $end; print OUT "\t"; print OUT $hash{$chr}{$start}{$end}[0]; print OUT "\t"; print OUT $hash{$chr}{$start}{$end}[1]; print OUT "\t"; print OUT $hash{$chr}{$start}{$end}[2]; print OUT "\n";
                }
           }
}


              
                  
                
