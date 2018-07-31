#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=Index123_2_pass_junction_counts_intron_retention_more_than_5_in_all_samples_sorted_list.txt; data_file_2=output_long_intron.gff; data_file_3=output_short_intron.gff; data_file_4=UNION_junc_coor_with_junction_ID_more_than_5_reads_overlap_from_all_samples.txt; data_file_5=For_DEXSeq_intron_retention_splicing_event.txt

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:preparing_intron_retention_event_for_DEXseq_combine.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4> <data_file_5>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4, $data_file_5)  = @ARGV;

my %hash_long; #$hash{junction_id}{left/right}=[chr,start,end];
my %hash_short;#$hash{junction_id}=[chr,strand,start,end];
my %hash_junc;#$hash{junction_id}=[chr,strand,start,end];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(IN3, "<$data_file_3") or die "can't open input2 file: $!";
open(IN4, "<$data_file_4") or die "can't open input2 file: $!";
open(OUT, ">$data_file_5") || die "can't open OUT file";

while(<IN2>) {
      chomp;
      my @long=split(/\s+/,$_);
      my @longid=split(/=/,$long[8]);
      my @longside=split(/-/,$longid[1]);
      $hash_long{$longside[0]}{$longside[1]}[0]=$long[0];
      $hash_long{$longside[0]}{$longside[1]}[1]=$long[3];
      $hash_long{$longside[0]}{$longside[1]}[2]=$long[4];
      #print $longside[0]; print "\t"; print $longside[1]; print "\n";
      }
close IN2;

while(<IN3>) {
      chomp;
      my @short=split(/\s+/,$_);
      my @shortid=split(/=/,$short[8]);
      $hash_short{$shortid[1]}[0]=$short[0];
      $hash_short{$shortid[1]}[1]=$short[3];
      $hash_short{$shortid[1]}[2]=$short[4];
      #print $shortid[1]; print "\n";
      }
close IN3;

while(<IN4>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash_junc{$array[4]}[0]=$array[0];
      $hash_junc{$array[4]}[1]=$array[3];
      $hash_junc{$array[4]}[2]=$array[1];
      $hash_junc{$array[4]}[3]=$array[2];
      #print $array[4]; print "\n";
}
close IN4;

while(<IN1>) {
      chomp;
      my $id=$_;
      my @junction=split(/-/,$id);
      my $left=1;
      my $right=2;
      #print $junction[0]; print "\t"; print $junction[1]; print "\n";
      if($junction[0] =~ /left/) {
          #print "left"; print "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[2]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[3];print OUT "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_long{$junction[1]}{$left}[1]; print OUT "\t"; print OUT $hash_long{$junction[1]}{$left}[2]; print OUT "\n";
      }
      if($junction[0] =~ /right/) {
          #print "right"; print "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[2]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[3];print OUT "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_long{$junction[1]}{$right}[1]; print OUT "\t"; print OUT $hash_long{$junction[1]}{$right}[2]; print OUT "\n"; 
     }
      if($junction[0] =~ /short/) {
          #print "short"; print "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[2]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[3];print OUT "\n";
          print OUT $id; print OUT "\t"; print OUT $hash_junc{$junction[1]}[0]; print OUT "\t"; print OUT $hash_junc{$junction[1]}[1]; print OUT "\t"; print OUT $hash_short{$junction[1]}[1]; print OUT "\t"; print OUT $hash_short{$junction[1]}[2]; print OUT"\n";
    } 

}

close IN1;



