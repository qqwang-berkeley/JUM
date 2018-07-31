#!/usr/bin/perl
use strict;
use warnings;

#data_file_1=INDEX3_fn_count.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:vary_for_replicate.pl <data_file_1> \n";

}

my ($data_file_1)  = @ARGV;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      #my $random_number = 0.05*(-1+rand(2));
      my $random_number = 0.1*(-1+rand(2));
      my $count = ($random_number+1)*$array[1];
      print $array[0]; print "\t"; print $array[1]; print "\t"; printf("%.0f", $count); print "\n";
}
      
           
