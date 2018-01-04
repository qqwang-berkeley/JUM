#!/usr/bin/perl
use strict;
use warnings;

#data_file_1=unique_AS_event_with_gene_names.txt; data_file_2=AS_XXXXXXX_sorted_with_dpsi; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:identify_gene_name_for_JUM_output_3.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") || die "can't open INPUT2 file";

my %hash;
my %final;
my @array;
my $origin;
my @temp;

while(<IN1>) {
      chomp;
      @array=split(/\s+/,$_);
      $hash{$array[0]}[0]=$array[1];
}

close IN1;

while(<IN2>) {
     chomp;
      $origin=$_;
      if($origin =~ /AS_event/) {
	      print "Gene"; print "\t"; print $origin; print "\n";
      }
      else {
           @temp=split(/\s+/,$_);
	   #if(exists $final{$temp[0]}) {
           #		   print " "; print "\t"; print $origin; print "\n";
	   #}
	   #else {
	   #	   $final{$temp[0]}[0]=0;
                   print $hash{$temp[0]}[0]; print "\t"; print $origin; print "\n";
          }
	  #}
}
close IN2;
              
                  
                
