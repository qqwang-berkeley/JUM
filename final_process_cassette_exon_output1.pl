#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
#data_file_1=cassette_sorted.txt; $data_file2=$condition_num; $data_file_3=$ctrl_num; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:final_process_cassette_exon_output.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my $count=0;
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") || die "can't open OUT1 file";
open(OUT2, ">$data_file_3") || die "can't open OUT2 file";

my @id=[0,0,0,0];
my $coor;
my @pvalue=[0,0,0,0];
my @qvalue=[0,0,0,0];
my @array;
my $pmin;
my $qmin;
while(<IN1>) {
      chomp;
      if($_ !~ /raw_count/) {
              if($count % 4 == 1) {
	      @array=split(/\s+/,$_);
	      $id[0]=$array[12];
	      $id[3]=$array[13];
	      $pvalue[0]=$array[6];
	      $qvalue[0]=$array[7];
	      $count=$count+1;
              } 
	      elsif($count % 4 == 2) {
              @array=split(/\s+/,$_);
	      $id[2]=$array[12];
	      $pvalue[1]=$array[6];
	      $qvalue[1]=$array[7];
	      $count=$count+1;
              }
              elsif($count % 4 == 3) {
              @array=split(/\s+/,$_);
	      $id[1]=$array[13]+1;
	      $pvalue[2]=$array[6];
	      $qvalue[2]=$array[7];
	      $count=$count+1;
              }
              elsif($count % 4 == 0) {
              @array=split(/\s+/,$_);
	      $pvalue[3]=$array[6];
	      $qvalue[3]=$array[7];
	      $count=$count+1;
	      $coor = $array[11] . "_" . $array[15] . "_" . $id[0] . "_" . $id[1] . "_" . $id[2] . "_" . $id[3];
	      $pmin = min @pvalue;
	      s/NA/1/g for @qvalue;
	      $qmin = min @qvalue;
	      print OUT1 $array[0]; print OUT1 "\t"; print OUT1 $coor; print OUT1 "\t"; print OUT1 $array[11]; print OUT1 "\t"; print OUT1 $array[15]; print OUT1 "\t"; print OUT1 $id[0]; print OUT1 "\t"; print OUT1 $id[1]; print OUT1 "\t"; print OUT1 $id[2]; print OUT1 "\t"; print OUT1 $id[3]; print OUT1 "\t" ; print OUT1 $pmin; print OUT1 "\t"; print OUT1 $qmin; print OUT1 "\t"; print OUT1 $array[-1]; print OUT1 "\n";
	      print OUT2 $array[1]; print OUT2 "\t"; print OUT2 $coor; print OUT2 "\n";
             }
     }
      else {  @array=split(/\s+/,$_);
	      print OUT1 "Gene"; print OUT1 "\t"; print OUT1 "AS_event_ID"; print OUT1 "\t"; print OUT1 "chromosome"; print OUT1 "\t"; print OUT1 "strand"; print OUT1 "\t"; print OUT1 "upstream_exon_end_coor"; print OUT1 "\t"; print OUT1 "cassette_exon_start_coor"; print OUT1 "\t"; print OUT1 "cassette_exon_end_coor"; print OUT1 "\t"; print OUT1 "downstream_exon_start_coor"; print OUT1 "\t" ; print OUT1 "pvalue"; print OUT1 "\t"; print OUT1 "qvalue"; print OUT1 "\t"; print OUT1 "deltaPSI_"; print OUT1 $array[8]; print OUT1 "-"; print OUT1 $array[9]; print OUT1 "\n";
	      $count=$count+1;
}
}
