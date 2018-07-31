#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=${files%Aligned.out_coverage.bed}"_"coverage_temp_long_intron_overlap_"$pvalue_padj"_"$cutoff".txt; data_file_2=${files%Aligned.out_coverage.bed}"_"temp_long_intron_retention_junction_coordinate_with_read_num_"$pvalue_padj"_"$cutoff".txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:count_intron_read_long_intron_retention_step2.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %section; #$section{chr}{start}{end}=[#num];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

my $num;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_); 
              for($num=0; $num<=$array[7]-$array[6]-1; $num++) {
		      if(($array[6]+$num >= $array[1]) && ($array[6]+$num+1 <= $array[2])) {
		      print OUT $array[5]; print OUT "\t"; print OUT $array[6]+$num; print OUT "\t"; print OUT $array[6]+$num+1; print OUT "\t"; print OUT $array[3]; print OUT "\t"; print OUT $array[4]; print OUT "\t"; print OUT $array[8]; print OUT "\n";
                      }
              }
}

close IN1;

