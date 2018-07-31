#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Descriptive;
#data_file_1=temp_long_intron_retention_junction_coordinate_with_read_num.txt; date_file_2=long_intron_retention_with_linear_fitting.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:determining_rightful_long_intron_retention_event.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %junction; #$junction{junction_id}=[#_coor1][#_coor2]....;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

my $indicator=1;
my $id;
my @y=();
my $outlier_indicator=0;
my $stat;
my $q1;
my $q3;
my $q4;
my $q0;
my $median;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      if($indicator == 1) {
	      $id = $array[4];
	      push @y, $array[5];
	      $indicator = $indicator + 1;
      }
      else {
	      if($array[4] eq $id) {
		      push @y, $array[5];
	      }
	      else {
                   $outlier_indicator=0;
                   $stat = Statistics::Descriptive::Full->new();
                   $stat->add_data(\@y);
                   $q1 = $stat->quantile(1);
                   $q3 = $stat->quantile(3);
                   $q0 = $stat->quantile(0);
                   $q4 = $stat->quantile(4);
                   $median = $stat->median();
                      if (($q0 < $median/4) || ($q1 < $q3/3) || ($q4 > $median*4)) {
                      $outlier_indicator=1;
                      }              
                      print OUT $id; print OUT "\t"; print OUT $outlier_indicator; print OUT "\t";
                      if ( grep( /^0$/, @y ) ) {       #this is to determine if there is any zero readout for read mapped to the intron. If yes, not valid intron retention and will be marked as 1, otherwise 0
                      print OUT "1"; print OUT "\n";
                      }
                      else { print OUT "0"; print OUT "\n";
                      }
		   $id = $array[4];
		   @y = ();
		   push @y, $array[5];
               }
	       $indicator = $indicator + 1;
       }
}

close IN1;

                  $outlier_indicator=0;
                   $stat = Statistics::Descriptive::Full->new();
                   $stat->add_data(\@y);
                   $q1 = $stat->quantile(1);
                   $q3 = $stat->quantile(3);
                   $q0 = $stat->quantile(0);
                   $q4 = $stat->quantile(4);
                   $median = $stat->median();
                      if (($q0 < $median/4) || ($q1 < $q3/3) || ($q4 > $median*4)) {
                      $outlier_indicator=1;
                      }              
                      print OUT $id; print OUT "\t"; print OUT $outlier_indicator; print OUT "\t";
                      if ( grep( /^0$/, @y ) ) {       #this is to determine if there is any zero readout for read mapped to the intron. If yes, not valid intron retention and will be marked as 1, otherwise 0
                      print OUT "1"; print OUT "\n";
                      }
                      else { print OUT "0"; print OUT "\n";
                      }
       
