#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=more_than_five_profiled_total_AS_event_junction_first_processing_for_DEXSeq_reference_building.txt; date_file_2=UNION_junc_coor_with_junction_ID_more_than_5_reads_overlap_from_all_samples.txt; #data_file_3=output; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_1.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %junction; #$junction{chr}{strand}{start}{end}=[junction_ID];
#my %pattern; #$pattern{number}=[AS_event1, AS_event2, ..., ]
#my %metric;  #$metric{number}=[#_of_existence_AS_1, #_of_existence_AS_2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") or die "can't open input3 file: $!";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $junction{$array[0]}{$array[3]}{$array[1]}{$array[2]}[0]=$array[4];
}
close IN2;

while(<IN1>) {
      chomp;
      my @line=split(/\s+/,$_);
      my $temp=$line[0];
      my @name=split(/:/,$line[0]);
      print OUT $name[0]; print OUT "\t"; print OUT $temp; print OUT "\t"; print OUT $line[1]; print OUT "\t"; print OUT $line[2]; print OUT "\t"; print OUT $line[3]; print OUT "\t"; print OUT $line[4]; print OUT "\t"; print OUT $junction{$line[1]}{$line[2]}{$line[3]}{$line[4]}[0]; print OUT "\n"; 
}

close IN1; 

                  
       
