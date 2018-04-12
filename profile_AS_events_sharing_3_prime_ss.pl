#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=UNION_junc_coor.txt; data_file_2=union_junc_coor_3_prime_ss_list_wiht_alternative_5_ss.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profile_AS_events_sharing_3_prime_ss.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %junction; #$junction{junction_id}=[chr][strand][start][end];
my %sscoor; #$5sscoor{chromsome}{strand}{coor}=[junctionID1, junctionID2, etc];
my $q;
my $p;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $sscoor{$array[0]}{$array[1]}{$array[2]}[0]=0;
}

close IN2;

while(<IN1>) {
      chomp;
      my @temp=split(/\s+/,$_);
      if (exists $sscoor{$temp[0]}{$temp[1]}{$temp[3]}) {
           $sscoor{$temp[0]}{$temp[1]}{$temp[3]}[0]=$sscoor{$temp[0]}{$temp[1]}{$temp[3]}[0]+1;
           $sscoor{$temp[0]}{$temp[1]}{$temp[3]}[$sscoor{$temp[0]}{$temp[1]}{$temp[3]}[0]]=$temp[4];
           }
      $junction{$temp[4]}[0]=$temp[0];
      $junction{$temp[4]}[1]=$temp[1];
      $junction{$temp[4]}[2]=$temp[2];
      $junction{$temp[4]}[3]=$temp[3];
}

close IN1;

foreach my $i (keys %sscoor) {
      foreach my $j (keys %{$sscoor{$i}}) {
             foreach my $k (keys %{$sscoor{$i}{$j}}) {
                       for($p=1; $p<=$sscoor{$i}{$j}{$k}[0]; $p++) {
                             print "3"; 
                             for($q=1; $q<=$sscoor{$i}{$j}{$k}[0]; $q++) {
                                 print "_"; print $sscoor{$i}{$j}{$k}[$q]; 
                             }
                             print "\t";
                             print $junction{$sscoor{$i}{$j}{$k}[$p]}[0]; print "\t"; print $junction{$sscoor{$i}{$j}{$k}[$p]}[1]; print "\t"; print $junction{$sscoor{$i}{$j}{$k}[$p]}[2]; print "\t"; print $junction{$sscoor{$i}{$j}{$k}[$p]}[3]; print "\n";
                       }
            }
      }
}                           

