#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=AS_differential_pvalue_0.05_A5SS_event_AS_structure_list.txt; data_file_2=AS_differential_JUM_output.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_4.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %AS_structure; #$AS_structure{AS_structure_id}=[0];
#my %AS_event; #$AS_event{AS_event_ID_linked_by_*}=[# of total AS structures, AS structure 1, AS structure 2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
#open(IN3, "<$data_file_3") or die "can't open input2 file: $!";
#open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

while(<IN1>) {
      chomp;
       $AS_structure{$_}[0]=0;
}
close IN1;

my $j;
my $ori;
while(<IN2>) {
     chomp;
     $ori=$_;
     my @pattern=split(/\s+/,$_);
     if($ori =~ /AS_structure_ID/) {
          print "AS_event_ID"; print "\t"; print $ori; print "\n";
     }
     else {
          if(exists $AS_structure{$pattern[0]}) {
             print $pattern[0]; print "\t"; print $ori; print "\n";
          }
     }
}
         
     #my $pattern=$_;
     #my $ori=$_;
     #my @pattern=split(/\*/,$_);
     #$AS_event{$ori}[0]=0;
     #   for($j=0; $j<=$#pattern; $j++) {
     #        if(exists $AS_event{$pattern[$j]}) {
                
        #  foreach my $i (keys %AS_event) {
        #       if($pattern =~ /$i/) {
        #             print $pattern; print "\t"; print $i; print "\t"; print $AS_event{$i}[0]; print "\n";
        #       }
        #  }
close IN2;

                  
       
