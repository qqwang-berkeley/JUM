#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=test_intron_delta_consistant_list.txt; data_file_2=UNION_junc_coor_with_junction_ID_more_than_*; input3=output #

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:screening_intron_retention_event_for_no_inner_junction.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %hash; #$hash{$junctionID}=[chr,strand,start,end];
my %hash_pos_junc; #$hash_pos_junc{$junctionID}=[1];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") or die "can't open input3 file: $!";

my @array;
my @temp;
my $chr;
my $strand;
my $start;
my $end;
my $indi;
my $ori;

while(<IN2>) {
      chomp;
      @array=split(/\s+/,$_);
      $hash{$array[4]}[0]=$array[0];
      $hash{$array[4]}[1]=$array[3];
      $hash{$array[4]}[2]=$array[1];
      $hash{$array[4]}[3]=$array[2];
}
close IN2;

while(<IN1>) {
      chomp;
      $ori=$_;
      @array=split(/\s+/,$_);
      @temp=split(/-/,$array[0]);
      if($ori !~ /right/)  {
        $chr=$hash{$temp[1]}[0];
        $strand=$hash{$temp[1]}[1];
        $start=$hash{$temp[1]}[2];
        $end=$hash{$temp[1]}[3];
        $indi = 0;
           foreach my $id (keys %hash) {
		   if(($hash{$id}[0] eq $chr) && ($hash{$id}[1] eq $strand)) {
			   if((($hash{$id}[2] == $start) && ($hash{$id}[3] < $end)) || (($hash{$id}[2] > $start) && ($hash{$id}[3] == $end))) { 
                                    $indi = 1;
				    last;
			    }
		    }
	    }
          if($indi == 0) {
	      $hash_pos_junc{$temp[1]}[0]=1; print OUT $ori; print OUT "\n";
          }
       }
     else {
	     if (exists $hash_pos_junc{$temp[1]}) {
		     print OUT $ori; print OUT "\n";
	     }
     }
}


