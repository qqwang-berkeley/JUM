#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=AS_differential_pvalue_0.05_valid_long_intron_retention_list.txt; #data_file_2=AS_differential_pvalue_0.05_valid_short_intron_retention_list.txt; data_file_3=AS_differential_JUM_output.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_6.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %AS_event;
my %AS_event_long; #$AS_event_long{AS_event_id}{left}{sub-junction_ID}=[line]
my %AS_event_short; #$AS_event_short{AS_event_id}{sub-junction_ID}=[line];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(IN3, "<$data_file_3") or die "can't open input2 file: $!";

while(<IN1>) {
      chomp;
      $AS_event{$_}[0]=0;
}
close IN1;

while(<IN2>) {
      chomp;
      $AS_event{$_}[0]=0;
}
close IN2;

my $ori;
my $tem;
while(<IN3>) {
     chomp;
     $ori=$_;
     my @pattern=split(/\s+/,$_);
     #if($ori =~ /AS_structure_ID/) {
     #     print "AS_event_ID"; print "\t"; print $ori; print "\n";
     #}
          if($pattern[0] =~ /left/) {
               $tem=$pattern[0];
               $pattern[0] =~ s/intronretentionleft-//; 
               $AS_event_long{$pattern[0]}{$tem}{$pattern[1]}[0]=$ori;
          }
          if($pattern[0] =~ /right/) {
               $tem=$pattern[0];
               $pattern[0] =~ s/intronretentionright-//;    
               $AS_event_long{$pattern[0]}{$tem}{$pattern[1]}[0]=$ori;
          }
          if($pattern[0] =~ /short/) {
               $AS_event_short{$pattern[0]}{$pattern[1]}[0]=$ori;
          }        
}
close IN3;

my $long_ID;
foreach my $key (keys %AS_event) {
        if($key =~ /short/) {
               foreach my $sh (keys %{$AS_event_short{$key}}) {
                    ($long_ID = $key) =~ s/short//; 
                    print $long_ID; print "\t"; print $AS_event_short{$key}{$sh}[0]; print "\n";
               }
        }
        else {
              foreach my $lo (sort keys %{$AS_event_long{$key}}) {
                     foreach my $su (sort keys %{$AS_event_long{$key}{$lo}}) {
                          $long_ID="intronretention" . '_' . $key; 
                          print $long_ID; print "\t"; print $AS_event_long{$key}{$lo}{$su}[0]; print "\n";
                     }
              }
       }
} 


                  
       
