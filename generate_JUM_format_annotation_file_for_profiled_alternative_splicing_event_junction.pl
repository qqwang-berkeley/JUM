#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=non_zero_profiled_total_alternative_splicing_event_junction_second_processing_for_DEXSeq_reference_building.txt; data_file_2=non_zero_profiled_total_alternative_splicing_event_junction_first_processing_for_DEXSeq_reference_building.txt; data_file_3=non_zero_profiled_total_alternative_splicing_event_junction_DEXSeq_annotation.gff; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:generate_annotation_file_for_non_zero_profiled_alternative_splicing_event_junction_for_DEXSeq.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my $o=0;
my %hash; #$hash{gene_name}{chr}{strand}=[start][end];
my @temp;
my @name;
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      if(exists $hash{$array[0]}{$array[1]}{$array[2]}) { }
      else {
             $hash{$array[0]}{$array[1]}{$array[2]}[0]=$array[3];
             $hash{$array[0]}{$array[1]}{$array[2]}[1]=$array[4];
      }
}

close IN1;

while(<IN2>) {
      chomp;
      @temp=split(/\s+/,$_);
      @name=split(/:/,$temp[0]);
      if((exists $hash{$name[0]}{$temp[1]}{$temp[2]}) && ($name[1] eq "001")) {
          print OUT $temp[1]; print OUT "\t"; print OUT "splice_junction_refGene_mm10_and_merged_from_five_samples.gtf"; print OUT "\t"; print OUT "aggregate_gene"; print OUT "\t"; print OUT $hash{$name[0]}{$temp[1]}{$temp[2]}[0]; print OUT "\t"; print OUT $hash{$name[0]}{$temp[1]}{$temp[2]}[1]; print OUT "\t"; print OUT "."; print OUT "\t"; if (($temp[2] ne "+") && ($temp[2] ne "-")) {print OUT "*";} else {print OUT $temp[2];}  print OUT "\t"; print OUT "."; print OUT "\t"; print OUT "gene_id"; print OUT " "; print OUT $name[0]; print OUT "\n";
          print OUT $temp[1]; print OUT "\t"; print OUT "splice_junction_refGene_mm10_and_merged_from_five_samples.gtf"; print OUT "\t"; print OUT "exonic_part"; print OUT "\t"; print OUT $temp[3]; print OUT "\t"; print OUT $temp[4]; print OUT "\t"; print OUT "."; print OUT "\t"; if (($temp[2] ne "+") && ($temp[2] ne "-")) {print OUT "*";} else {print OUT $temp[2];} print OUT "\t"; print OUT "."; print OUT "\t"; print OUT "exonic_part_number"; print OUT " "; print OUT $name[1]; print OUT ";"; print OUT " "; print OUT "gene_id"; print OUT " "; print OUT $name[0]; print OUT "\n";  
      }

      if((exists $hash{$name[0]}{$temp[1]}{$temp[2]}) && ($name[1] ne "001")) {
          print OUT $temp[1]; print OUT "\t"; print OUT "splice_junction_refGene_mm10_and_merged_from_five_samples.gtf"; print OUT "\t"; print OUT "exonic_part"; print OUT "\t"; print OUT $temp[3]; print OUT "\t"; print OUT $temp[4]; print OUT "\t"; print OUT "."; print OUT "\t"; if (($temp[2] ne "+") && ($temp[2] ne "-")) {print OUT "*";} else {print OUT $temp[2];}  print OUT "\t"; print OUT "."; print OUT "\t"; print OUT "exonic_part_number"; print OUT " "; print OUT $name[1]; print OUT ";"; print OUT " "; print OUT "gene_id"; print OUT " "; print OUT $name[0]; print OUT "\n";

       }
        
}
       
               
close IN2;




