#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=count_table_index1_short_intron (generated from R); data_file_2=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples.txt; data_file_3=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_short_intron_retention.txt; data_file_4=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_short_intron_retention_list.txt; threshold # 

#more than thresholded # of reads

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:making_intron_retention_count_table_short_intron.pl <data_file_1> <data_file_2> <data_file_3><data_file_4> <threshold>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4, $threshold)  = @ARGV;

my $o=0;
my %hash; #$hash{junction_id}=[chr,strand,start,end,count_no_intron,count_with_intron]
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";
open(OUT2, ">$data_file_4") || die "can't open OUT2 file";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[4]}[0]=$array[0];
      $hash{$array[4]}[1]=$array[1];
      $hash{$array[4]}[2]=$array[2];
      $hash{$array[4]}[3]=$array[3];
      $hash{$array[4]}[4]=$array[5];
      $hash{$array[4]}[5]=0;
}

close IN2;

while(<IN1>) {
      chomp;
      my @temp=split(/\s+/,$_);
      if(exists $hash{$temp[0]}) {
                $hash{$temp[0]}[5]=$temp[1];
          }
}
close IN1;

foreach my $n (keys %hash) {
	#if($hash{$n}[5] >= $threshold) {
         print OUT "intronretentionshort-"; print OUT $n; print OUT ":"; print OUT "001"; print OUT "\t"; 
         print OUT $hash{$n}[4]; print OUT "\n";
         print OUT "intronretentionshort-"; print OUT $n; print OUT ":"; print OUT "002"; print OUT "\t";
         print OUT $hash{$n}[5]; print OUT "\n";
         if($hash{$n}[5] >= $threshold) {
	 print OUT2 "intronretentionshort-"; print OUT2 $n; print OUT2 "\n";
         }
}


