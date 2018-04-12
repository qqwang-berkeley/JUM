#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=non_zero_profiled_total_alternative_splicing_event_junction.txt; data_file_2=non_zero_profiled_total_alternative_splicing_event_junction_second_processing_for_DEXSeq_reference_building.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:second_processing_generate_DEXSeq_reference_format_for_non_zero_profiled_total_alternative_splicing_event_junction.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %hash; #$hash{junction_id}{chr}{strand}=[start,end];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") || die "can't open OUT file";
while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      if(exists $hash{$array[0]}{$array[1]}{$array[2]}) {
            if($array[3] < $hash{$array[0]}{$array[1]}{$array[2]}[0]) {
                 $hash{$array[0]}{$array[1]}{$array[2]}[0]=$array[3];
            }
            if($array[4] > $hash{$array[0]}{$array[1]}{$array[2]}[1]) {
                 $hash{$array[0]}{$array[1]}{$array[2]}[1]=$array[4];
            }
      }
      else {
         $hash{$array[0]}{$array[1]}{$array[2]}[0]=$array[3];
         $hash{$array[0]}{$array[1]}{$array[2]}[1]=$array[4];
      }
}

close IN1;

 foreach my $id (keys %hash) {
             foreach my $chr (keys %{$hash{$id}}) {
                   foreach my $strand (keys %{$hash{$id}{$chr}}) {
                         print OUT $id; print OUT "\t"; print OUT $chr; print OUT "\t"; print OUT $strand; print OUT "\t"; print OUT $hash{$id}{$chr}{$strand}[0]; print OUT "\t"; print OUT $hash{$id}{$chr}{$strand}[1]; print OUT "\n"
                    }
             }
}
