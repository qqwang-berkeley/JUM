#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=non_zero_profiled_total_alternative_splicing_event_junction_first_processing_for_JUM_reference_building.txt; data_file_2=ctrl1_remap_junction_counts.txt; data_file_3=CTRL1_fn_count.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:making_junction_count_file_matching_annotation_for_JUM_in_intron_retention.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my $o=0;
my %hash; #$hash{chr}{strand}{start}{end}=[count_number];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[0]}[0]=$array[1];
}

close IN2;


while(<IN1>) {
      chomp;
      my @temp=split(/\s+/,$_);
      if(exists $hash{$temp[0]}) {
           print OUT $temp[0]; print OUT "\t"; print OUT $hash{$temp[0]}[0]; print OUT "\n";
      }
      else {
           print OUT $temp[0]; print OUT "\t"; print OUT $o; print OUT "\n";
      }
}
close IN1;

       



