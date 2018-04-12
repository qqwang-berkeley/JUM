#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=non_zero_profiled_total_alternative_splicing_event_junction.txt; data_file_2=non_zero_profiled_total_alternative_splicing_event_junction_first_processing_for_DEXSeq_reference_building.txt; 
#this script is written that requires the input file to be sorted by gene name first
my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:first_processing_generate_DEXSeq_reference_format_for_non_zero_profiled_total_alternative_splicing_event_junction.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my $o=1;
my $number;
my $con;
my %hash; #$hash{gene_name}{chr}{strand}=[$o];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") || die "can't open OUT file";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      if(exists $hash{$array[0]}{$array[1]}{$array[2]}) {
             $hash{$array[0]}{$array[1]}{$array[2]}[0]=$hash{$array[0]}{$array[1]}{$array[2]}[0]+1;
             print OUT $array[0]; print OUT ":"; printf OUT '%03s',$hash{$array[0]}{$array[1]}{$array[2]}[0]; print OUT "\t"; print OUT $array[1]; print OUT "\t"; print OUT $array[2]; print OUT "\t"; print OUT $array[3]; print OUT "\t";print OUT $array[4]; print OUT "\n";
      }
      else { 
             $hash{$array[0]}{$array[1]}{$array[2]}[0]=$o;
             print OUT $array[0]; print OUT ":";  printf OUT '%03s',$hash{$array[0]}{$array[1]}{$array[2]}[0]; print OUT "\t"; print OUT $array[1]; print OUT "\t"; print OUT $array[2]; print OUT "\t"; print OUT $array[3]; print OUT "\t";print OUT $array[4]; print OUT "\n"; 
      }
}

close IN1;




