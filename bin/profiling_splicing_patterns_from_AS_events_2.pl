#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=temp_long_intron_retention_junction_coordinate_with_read_num.txt; date_file_2=long_intron_retention_with_linear_fitting.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_2.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %junction; #$junction{junction_id}=[# of existence in all AS events];
#my %pattern; #$pattern{number}=[AS_event1, AS_event2, ..., ]
#my %metric;  #$metric{number}=[#_of_existence_AS_1, #_of_existence_AS_2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
       if(exists $junction{$array[6]}) {
          $junction{$array[6]}[0]=$junction{$array[6]}[0]+1;
          $junction{$array[6]}[$junction{$array[6]}[0]]=$array[0];
       }
       else {
          $junction{$array[6]}[0]=1;
          $junction{$array[6]}[1]=$array[0];
       }
}
close IN1;

my $l;
foreach my $i (keys %junction) {
       print OUT $i; print OUT "\t"; print OUT $junction{$i}[0]; 
          for($l=1; $l<=$junction{$i}[0]; $l++) {
                 print OUT "\t"; print OUT $junction{$i}[$l]; 
          }
       print OUT "\n";
}
                  
       
