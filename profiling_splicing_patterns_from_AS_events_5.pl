#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=total_cassette_exon_event.txt; data_file_2=AS_differential_JUM_output.txt;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events_5.pl <data_file_1> <data_file_2>\n";

}

my ($data_file_1, $data_file_2)  = @ARGV;

my %AS_structure; #$AS_structure{AS_structure_id}{$sub_junction_id}=[whole line];
my %AS_event; #$AS_event{AS_event_ID_linked_by_*}=[AS structure 1, AS structure 2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
#open(IN3, "<$data_file_3") or die "can't open input2 file: $!";
#open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

my $i;
while(<IN1>) {
      chomp;
      my $ori=$_;
      my @pattern=split(/\*/,$ori);
       for($i=0; $i<=$#pattern; $i++) {
            $AS_event{$ori}[$i]=$pattern[$i];
       }
}
close IN1;   
      #print $pattern[0]; print "\n"; print $pattern[1]; print "\n";
      #my @array=split(/\s+/,$_);
       #if(exists $junction{$array[6]}) {
       #   $junction{$array[6]}[0]=$junction{$array[6]}[0]+1;
       #   $junction{$array[6]}[$junction{$array[6]}[0]]=$array[0];
       #}
      #$AS_event{$array[0]}[0]=$array[1];
       #$AS_structure{$_}[0]=0;

while(<IN2>) {
     chomp;
     my $line=$_;
     my @string=split(/\s+/,$line);
     $AS_structure{$string[0]}{$string[1]}[0]=$line;
}

close IN2;

my $j;
foreach my $key (keys %AS_event) {
     for($j=0; $j<=$#{ $AS_event{$key}}; $j++) {
                 if(exists $AS_structure{$AS_event{$key}[$j]}) {
                     foreach my $m (keys %{$AS_structure{$AS_event{$key}[$j]}}) {
                        print $key; print "\t"; print $AS_structure{$AS_event{$key}[$j]}{$m}[0]; print "\n";
                     } 
                 }
     }
}

#if($pvalue_adj =~ /padj/) {
#    while(<IN2>) {
#      chomp;
#      if($_ =~ /AS_structure_ID/) { print "AS_event_ID"; print "\t"; print $_; print "\n";}
#      else {
#       my $cha=$_;
#       my @string=split(/\s+/,$cha);
#           if($string[6] <= $value) {
#                 foreach my $j (keys %AS_event) {
#                     if (grep( /$string[0]/, @{ $AS_event{$j}})) { 
#                          print $j; print "\t"; print $cha; print "\n";
#                     }
#                 }
#           }
#     }
#   }
#   close IN2;                 
#}  

