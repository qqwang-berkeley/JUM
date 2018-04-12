#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
#use Array::Compare;
#data_file_1=tt; date_file_2=t1; #data_file_3=output_A5SS; $data_file_4=output_A3SS; $data_file_5=output_cassette_exon; $data_file_6=output_complicated; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_3.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4> <data_file_5> <data_file_6>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4, $data_file_5, $data_file_6)  = @ARGV;

my %junction; #$junction{ID}=[occurance];
my %AS_component; #$AS_component{AS_ID}=[#ofsubASevents, sub_event1, sub_event2, ..., ] 
my %AS_occurance; #$AS_occurance{AS_ID}=[#ofsubASevents,occurance_sub_event1, occurance_sub_event2, ..., ]
my %AS_metric; #$AS_metric{AS_ID}=[number of 2s; number of 1s];
my %AS_relations; #$AS_relations{AS_ID}=[related_subAS_1, related_subAS_2, ..., ]
my %pattern; #$pattern{ID}=[#AS_1, AS_2, ... , ]
my %mark; #$mark{AS_ID}=[1]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT1, ">$data_file_3") or die "can't open input3 file: $!";
open(OUT2, ">$data_file_4") or die "can't open input4 file: $!";
open(OUT3, ">$data_file_5") or die "can't open input5 file: $!";
open(OUT4, ">$data_file_6") or die "can't open input6 file: $!";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $junction{$array[0]}[0]=$array[1];
         if($array[1] == 1) {
            if(exists $AS_relations{$array[2]}) {
                 push @{ $AS_relations{$array[2]} }, $array[2];
            }
            else {
                 $AS_relations{$array[2]}[0] = $array[2];
            }
         }
         if($array[1] == 2) {
            if(exists $AS_relations{$array[2]}) {
                 push @{ $AS_relations{$array[2]} }, $array[3];
           }
            else {
                 $AS_relations{$array[2]}[0] = $array[3];
          }
            if(exists $AS_relations{$array[3]}) {
                 push @{ $AS_relations{$array[3]} }, $array[2];
           }
            else {
                 $AS_relations{$array[3]}[0] = $array[2];
          }
        }
}
close IN2;

#foreach my $ID (keys %AS_relations) {
#        my @temp = uniq @{ $AS_relations{$ID} };
#        @{ $AS_relations{$ID} } = @temp;
#        print OUT3 $ID; print OUT3 "\t"; 
#        foreach (@{ $AS_relations{$ID} }) {
#              print OUT3"$_\t"; 
#        }
#     print OUT3 "\n"; 
#} 


while(<IN1>) {
      chomp;
      my @line=split(/\s+/,$_);
          if(exists $AS_component{$line[0]}) {
              $AS_component{$line[0]}[0]=$line[3];
              $AS_occurance{$line[0]}[0]=$AS_occurance{$line[0]}[0]+1;
              $AS_component{$line[0]}[$AS_occurance{$line[0]}[0]] = $line[6];
              $AS_occurance{$line[0]}[$AS_occurance{$line[0]}[0]] = $junction{$line[6]}[0];                            

              #push @{ $AS_component{$line[0]} }, $line[6];
              #push @{ $AS_occurance{$line[0]} }, $junction{$line[6]}[0];
          }
          else {
              $AS_component{$line[0]}[0]=$line[3];
              $AS_occurance{$line[0]}[0]=1;
              $AS_component{$line[0]}[1] = $line[6];
              $AS_occurance{$line[0]}[1] = $junction{$line[6]}[0];
          } 
}

close IN1; 
      
foreach my $id (keys %AS_occurance) {
        my @AS_1 = @{ $AS_occurance{$id} }[1 .. $AS_occurance{$id}[0]];
        $AS_metric{$id}[0] = grep (/^2$/, @AS_1);
        $AS_metric{$id}[1] = grep (/^1$/, @AS_1);
        #$AS_metric{$id}[0] = grep (/^2$/, @{ $AS_occurance{$id} });       
        #$AS_metric{$id}[1] = grep (/^1$/, @{ $AS_occurance{$id} });
}

my $t;
my @temp_ori;
my @temp;
my @temp_after;
my $m;
my @temp_mid;

foreach my $pin (keys %AS_relations) {
         @temp_ori = sort @{ $AS_relations{$pin} };
         @temp_mid=@temp_ori;
             for($t=0; $t<=$#{ $AS_relations{$pin} }; $t++) {
                 push (@temp_ori, @{ $AS_relations{$AS_relations{$pin}[$t]} });
             }
         @temp = uniq @temp_ori;
         @temp_after = sort @temp;
             while(!(@temp_mid ~~ @temp_after)) {
                 @temp_mid = @temp_after;
                    for($t=0; $t<=$#temp_mid; $t++) {
                        push (@temp_after, @{ $AS_relations{$temp_mid[$t]} });
                    }
                @temp = uniq @temp_after;
                @temp_after = sort @temp;
            }
            my $PIN = join('*', @temp_after);
                if(exists $pattern{$PIN}) { }
                else {
                     @{ $pattern{$PIN} } = @temp_after;
                     if($PIN !~ /\*/) {
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^5_/) && ($AS_component{$PIN}[0] eq "+")) {
                                  print OUT2 $PIN; #print OUT2 "\t"; foreach (@temp_after) { print OUT2 "$_\t";} 
                                  print OUT2 "\n";
                              }
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^3_/) && ($AS_component{$PIN}[0] eq "-")) {
                                  print OUT2 $PIN; #print OUT2 "\t"; foreach (@temp_after) { print OUT2 "$_\t";} 
                                  print OUT2 "\n";
                              }
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^5_/) && ($AS_component{$PIN}[0] eq "0")) {
                                  print OUT2 $PIN; #print OUT2 "\t"; foreach (@temp_after) { print OUT2 "$_\t";} 
                                  print OUT2 "\n";
                              }
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^5_/) && ($AS_component{$PIN}[0] eq "-")) {
                                  print OUT1 $PIN; #print OUT1 "\t"; foreach (@temp_after) { print OUT1 "$_\t";} 
                                  print OUT1 "\n";
                              }
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^3_/) && ($AS_component{$PIN}[0] eq "+")) {
                                  print OUT1 $PIN; #print OUT1 "\t"; foreach (@temp_after) { print OUT1 "$_\t";} 
                                  print OUT1 "\n";
                              }
                             if(($AS_metric{$PIN}[0] == 0) && ($PIN =~ /^3_/) && ($AS_component{$PIN}[0] eq "0")) {
                                  print OUT1 $PIN; #print OUT1 "\t"; foreach (@temp_after) { print OUT1 "$_\t";} 
                                  print OUT1 "\n";
                              }
                   }
                    else {  my $index_zero_1_AS=0;
                            my $index_only_1_other_2_AS=0;
                            my $p;
                              for($p=0; $p<=$#temp_after; $p++) {
                                     if(($AS_metric{$temp_after[$p]}[0] >= 1) && ($AS_metric{$temp_after[$p]}[1] == 0)) {
                                               $index_zero_1_AS=$index_zero_1_AS+1;
                                     }
                                     if(($AS_metric{$temp_after[$p]}[0] >= 1) && ($AS_metric{$temp_after[$p]}[1] == 1)) {
                                               $index_only_1_other_2_AS=$index_only_1_other_2_AS+1;
                                    }
                              }
			      #if(($index_only_1_other_2_AS == 2) && ($#temp_after+1 == $index_zero_1_AS+$index_only_1_other_2_AS) && ($index_zero_1_AS > 0)) {
			      # if(($index_only_1_other_2_AS == 2) && ($#temp_after+1 == $index_zero_1_AS+$index_only_1_other_2_AS) && ($index_zero_1_AS > 1)) {
			      #	   print OUT4 $PIN;
			      #     print OUT4 "\n";
                              # }
			      #elsif(($index_only_1_other_2_AS == 2) && ($#temp_after+1 == $index_zero_1_AS+$index_only_1_other_2_AS) && ($index_zero_1_AS == 0)) {
                           if(($index_only_1_other_2_AS == 2) && ($#temp_after+1 == $index_zero_1_AS+$index_only_1_other_2_AS) && ($index_zero_1_AS == 0)) {

                           #if(($#temp_after == 1) && ($AS_occurance{$temp_after[0]}[0] == 2) && ($AS_occurance{$temp_after[1]}[0] == 2) && ($AS_metric{$temp_after[0]}[0] == 1) && ($AS_metric{$temp_after[0]}[1] == 1) && ($AS_metric{$temp_after[1]}[0] == 1) && ($AS_metric{$temp_after[1]}[1] == 1)) {
                                 print OUT3 $PIN;
                                 print OUT3 "\n";
                           }
                           else {
                                 print OUT4 $PIN;
                                 print OUT4 "\n";
                           }
                       } 
               }
               
}

#foreach my $try (keys %pattern) {
#           print $try; print "\t";
#           foreach (@{ $pattern{$try} }) {
#              print "$_\t";
#           }
#     print "\n";
#}
                  
       
