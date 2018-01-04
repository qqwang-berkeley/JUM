#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
#data_file_1=IR_short_sorted.txt; $data_file2=$condition_num; $data_file_3=$ctrl_num; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:output_reformat_2.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $condition_num, $ctrl_num)  = @ARGV;

my $o=0;
my %hash; #$hash{chr}{strand}{start}{end}=[count_number];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";

my @enter_1;
$#enter_1=$ctrl_num+$condition_num;
my @enter_2;
$#enter_2=$ctrl_num+$condition_num;

my @ratio;
$#ratio=$ctrl_num+$condition_num;

my $count=0;
my $indicator_1;
my $indicator_2;
my $dpsi;
my $av1;
my $av2;
my $sum1;
my $sum2;
my $s;
my $m;

while(<IN1>) {
      chomp;
    if($_ !~ /raw_count/) {
      $count=$count+1;
      my $ori=$_;
      my @array=split(/\s+/,$_);
         if($count % 2 == 1) {
		 for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
                 $enter_1[$s]=$array[15+$s];
	        }
	   print $ori; print "\n";
         } 
         if($count % 2 == 0) {
                 for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
                 $enter_2[$s]=$array[15+$s];
                 }
	 $indicator_1=0;
	 $indicator_2=0;
	 $sum1=0;
	 $sum2=0;
         for($s=0;$s<=$condition_num-1;$s++) {
	     if($enter_1[$s]+$enter_2[$s] != 0) {	
		    $ratio[$s]=$enter_2[$s]/($enter_1[$s]+$enter_2[$s]);
                    $indicator_1=$indicator_1+1;
	    }
	     else {
		     $ratio[$s]=0;
	     }
          }
          for($m=$s;$m<=$condition_num+$ctrl_num-1;$m++) {
	     if($enter_1[$m]+$enter_2[$m] != 0) {	
		    $ratio[$m]=$enter_2[$m]/($enter_1[$m]+$enter_2[$m]);
                    $indicator_2=$indicator_2+1;
	    }
	     else {
		     $ratio[$m]=0;
	     }
          }
         if(($indicator_1 != 0) && ($indicator_2 != 0)) {
                 for($s=0;$s<=$condition_num-1;$s++) {
			 $sum1=$sum1+$ratio[$s];
		 }
		 $av1=$sum1/$indicator_1;
		 for($m=$s;$m<=$condition_num+$ctrl_num-1;$m++) {
			 $sum2=$sum2+$ratio[$m];
		 }
		 $av2=$sum2/$indicator_2;
		 $dpsi=$av1-$av2;
	 }
	 else {
		 $dpsi="INF";
	 }
	 print $ori; print "\t"; print $dpsi; print "\n";
     }
     }
     else {
	     print $_; print "\t"; print "deltaPSI"; print "\n";
     }
}

