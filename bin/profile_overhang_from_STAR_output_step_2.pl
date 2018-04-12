#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=DRQW1A_2nd_pass_junction_counts.txt; data_file_2=Profiled_Index1_Aligned.out.spanning_junction_reads_junction_and_overhang_and_mapped_number.txt; data_file_3=output; data_file_4=back_output;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profile_overhang_from_STAR_output_step_2_no_strand_considered.pl <data_file_1> <data_file_2> <data_file_3><data_file_4>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4)  = @ARGV;

my $o=0;
my %hash; #$hash{chr}{start}{end}=[mapped_number,leftoverhang,rightoverhang];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";
open(OUT_1, ">$data_file_4") || die "can't open OUT file";

while(<IN2>) {
      chomp;
      my @junction=split(/\s+/,$_);
      $hash{$junction[0]}{$junction[1]}{$junction[2]}[0]=$junction[3];
      $hash{$junction[0]}{$junction[1]}{$junction[2]}[1]=$junction[4];
      $hash{$junction[0]}{$junction[1]}{$junction[2]}[2]=$junction[5];
}

while(<IN1>) {
      chomp;
      my $temp=$_;
      my @array=split(/\s+/,$_);
      $array[2]=$array[2]+1;
      $array[3]=$array[3]+1;
      if(exists $hash{$array[0]}{$array[2]}{$array[3]}) {
           print OUT $temp; print OUT "\t"; print OUT $hash{$array[0]}{$array[2]}{$array[3]}[0]; print OUT "\t"; print OUT $hash{$array[0]}{$array[2]}{$array[3]}[1]; print OUT "\t"; print OUT $hash{$array[0]}{$array[2]}{$array[3]}[2]; print OUT "\n";
      }
      else {
           print OUT_1 $temp; print OUT_1 "\n";
     }
}

              
                  
                
