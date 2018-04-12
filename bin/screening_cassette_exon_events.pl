#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);
#data_file_1=total_cassette_exon_event.txt ; date_file_2=reconstruct_splicing_pattern_input_1.txt ; data_file_3=cassette_exon.bed; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:screening_cassette_exon_events_1.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %junction; #$junction{juntion_ID}=[chr,strand,start,end];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") or die "can't open input3 file: $!";

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      $junction{$array[6]}[0]=$array[2];
      $junction{$array[6]}[1]=$array[3];
      $junction{$array[6]}[2]=$array[4];
      $junction{$array[6]}[3]=$array[5];
}
close IN2;

my @cassette;
my @A5SS;
my @A3SS;
my @common;
while(<IN1>) {
      chomp;
      my $ori=$_;
      my @cassette=split(/\*/,$_);
      my @SS3=split(/_Junction_/,$cassette[0]);
      my @SS5=split(/_Junction_/,$cassette[1]);
      my $A3SS1="Junction_" . $SS3[1];
      my $A3SS2="Junction_" . $SS3[2];
      my $A5SS1="Junction_" . $SS5[1];
      my $A5SS2="Junction_" . $SS5[2];
      #print $A3SS1; print "\n"; print $A3SS2; print "\n"; print $A5SS1; print "\n"; print $A5SS2; print "\n";
      @A3SS=qw();
      push @A3SS, $A3SS1; push @A3SS, $A3SS2;
      @A5SS=qw();
      push @A5SS, $A5SS1; push @A5SS, $A5SS2;
      #foreach (@A3SS) { print "$_\t";} print "\n"; foreach (@A5SS) { print "$_\t";} print "\n";
      @common= intersect(@A5SS, @A3SS);
      #print $common[0]; print "\n";
      #print $cassette[0]; print "\n";
      $cassette[0]=~s/^3_//;
      $cassette[0]=~s/_$common[0]//;
      $cassette[0]=~s/$common[0]_//;
      #$cassette[0]=~s/$common[0]//;
      #print $cassette[0]; print "\n";
      #print $cassette[1]; print "\n";
      $cassette[1]=~s/^5_//;
      $cassette[1]=~s/_$common[0]//; 
      $cassette[1]=~s/$common[0]_//;
      #$cassette[1]=~s/$common[0]//;
      #print $cassette[1]; print "\n";
      #if($junction{$cassette[1]}[3] <= $junction{$cassette[0]}[2]) {
      #if($junction{$cassette[1]}[3] < $junction{$cassette[0]}[2]) {
      if($junction{$cassette[1]}[3] < $junction{$cassette[0]}[2]-1) {

	      #print OUT $junction{$cassette[0]}[0]; print OUT "\t"; print OUT "cassette_exon"; print OUT "\t"; print OUT "cassette"; print OUT "\t"; print OUT $junction{$cassette[1]}[3]; print OUT "\t"; print OUT $junction{$cassette[0]}[2]; print OUT "\t"; print OUT "."; print OUT "\t"; print OUT $junction{$cassette[0]}[1]; print OUT "\t"; print OUT "."; print OUT "\t"; print OUT "cassette_id="; print OUT $_; print OUT "\n";             
              print OUT $junction{$cassette[0]}[0]; print OUT "\t"; print OUT $junction{$cassette[1]}[3]+1; print OUT "\t"; print OUT $junction{$cassette[0]}[2]; print OUT "\t"; print OUT $junction{$cassette[0]}[1]; print OUT "\t"; print OUT $_; print OUT "\n"; 
      }
}

close IN1; 

                  
       
