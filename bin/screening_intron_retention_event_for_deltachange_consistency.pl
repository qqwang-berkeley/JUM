#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
#data_file_1=non_zero_profiled_total_alternative_splicing_event_junction_first_processing_for_DEXSeq_reference_building.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:screening_intron_retention_event_for_deltachange_consistency.pl <data_file_1>\n";

}

my ($data_file_1)  = @ARGV;

my $o=0;
my %hash; #$hash{chr}{strand}{start}{end}=[count_number];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";

my @enter;
$#enter=4;
my @array;
my $left_sign1=0;
my $left_sign2=0;
my $right_sign1=0;
my $right_sign2=0;
my $count=0;

while(<IN1>) {
      chomp;
      if($_ !~ /AS_structure/) {
	   if($_ !~ /short/) { 
                $count=$count+1;
                @array=split(/\s+/,$_);
                    if($count % 4 == 1) {
			 @enter=();
			 $enter[0]=$array[0];
			 $left_sign1=0;
			 $left_sign2=0;
			 $right_sign1=0;
			 $right_sign2=0;
		           if($array[2] > 0) {
		              $left_sign1=1;
			   }
			   else {
				   $left_sign1=-1;
			   }
		   }
                   if($count % 4 == 2) {
                          $enter[1]=$array[0];
			   if($array[2] > 0) {
			      $left_sign2=1;
			   }
			   else {
				  $left_sign2=-1;
			  }
		  }
                  if($count % 4 == 3) {
                          $enter[2]=$array[0];
		           if($array[2] > 0) {
			      $right_sign1=1;
		          }
		           else {
			         $right_sign1=-1;
		         }
	         }
                  if($count % 4 == 0) {
                          $enter[3]=$array[0];
		           if($array[2] > 0) {
				  $right_sign2=1;
			   }
			   else {
				  $right_sign2=-1;
			   }
			   if((($left_sign1 == 1) && ($right_sign1 == 1) && ($left_sign2 == -1) && ($right_sign2 == -1)) || (($left_sign1 == -1) && ($right_sign1 == -1) && ($left_sign2 == 1) && ($right_sign2 == 1))) {
				  print $enter[0]; print "\n"; print $enter[1]; print "\n"; print $enter[2]; print "\n"; print $enter[3]; print "\n";
		          }
		  }
	  }
	  else {  @array=split(/\s+/,$_);
		  print $array[0]; print "\n";
	  }
  }
}



   
   
   
   
       



