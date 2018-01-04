#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=Index1_SJ.out.tab_strand_symbol_scaled; data_file_2=UNION_junc_coor.txt; data_file_2=Index1_second_pass_idxstats_junction.txt;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:prepare_count_file_after_STAR_2_pass.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %hash; #$hash{chr}{strand}{coor1}{coor2}=[junctionID];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input1 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";
my $i;

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
         $hash{$array[0]}{$array[3]}{$array[1]}{$array[2]}[0]=$array[4];
}

while(<IN1>) {
      chomp;
      my $line=$_;
      my @coor=split(/\s+/,$_);
           if(exists $hash{$coor[0]}{$coor[3]}{$coor[1]}{$coor[2]}) { 
               print OUT $line; print OUT "\t"; print OUT $hash{$coor[0]}{$coor[3]}{$coor[1]}{$coor[2]}[0]; print OUT "\n";
           }
           else {
              print OUT $line; print OUT "\t"; print OUT "NONE"; print OUT "\n";
          }
}


              
                  
                
