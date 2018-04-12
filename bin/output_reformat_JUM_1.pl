#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
#data_file_1=A5SS/A3SS/MXE/mixed_sorted.txt; data_file_2=condition_num; data_file_3=ctrl_num; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:output_reformat_JUM_1.pl <data_file_1> \n";

}

my ($data_file_1, $condition_num, $ctrl_num)  = @ARGV;

my $fold;
my $c;
my $av1;
my $av2;
my $sum1;
my $sum2;
my $s1;
my $s2;
my $dpsi;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";

while(<IN1>) {
      chomp;
      if($_ !~ /raw_count/) {
      my $ori=$_;
      my @array=split(/\s+/,$_);
      for($c=0; $c<$ctrl_num+$condition_num; $c++) {
           $array[15+$ctrl_num+$condition_num+$c] =~ s/%//;
      }
      $sum1=0;
      $sum2=0;
      for($s1=0; $s1<=$condition_num-1;$s1++) {
	      $sum1=$sum1+$array[15+$ctrl_num+$condition_num+$s1];
      }
      $av1=$sum1/$condition_num/100;
      for($s2=$s1; $s2<$ctrl_num+$condition_num; $s2++) {
	      $sum2=$sum2+$array[15+$ctrl_num+$condition_num+$s2];
      }
      $av2=$sum2/$ctrl_num/100;
      $dpsi=$av1-$av2;
      print $ori; print "\t"; print $dpsi; print "\n";
      }
      else {
        print $_; print "\t"; print "deltaPSI"; print "\n";
      }

}

close IN1;


       



