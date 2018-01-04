#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );

#data_file_1=refFlat.txt; data_file_2=AS_X_sorted_with_dpsi.txt; data_file_3=output_with_gene_names; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:identify_gene_name_for_JUM_output.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %hash; #$hash{composite_id}=[gene, transcript, chr, strand, start,end];
#my %AS; #$AS{AS_event_id}=[gene1,gene2,gene3,...];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";

my @array;
my $string;

while(<IN1>) {
      chomp;
      @array=split(/\s+/,$_);
      $string=join('_', $array[0], $array[1], $array[2], $array[3], $array[4]);
      $hash{$string}[0]=$array[0];
      $hash{$string}[1]=$array[1];
      $hash{$string}[2]=$array[2];
      $hash{$string}[3]=$array[3];
      $hash{$string}[4]=$array[4];
      $hash{$string}[5]=$array[5];
}
close IN1;

my $tagg;
my @name;

while(<IN2>) {
   chomp;
   my $origin=$_;
   my @junctionID=split(/\s+/,$_);
   $tagg=0;
   @name=();
   if($origin !~ /AS_event/) {
      foreach my $id (keys %hash) {
                        if(($junctionID[10] eq $hash{$id}[2]) && ($junctionID[14] eq $hash{$id}[3])) {
                             if ((($junctionID[11] <= $hash{$id}[4]) && ($junctionID[12] >= $hash{$id}[4]) && ($junctionID[12] <= $hash{$id}[5])) || (($junctionID[11] >= $hash{$id}[4]) && ($junctionID[11] <= $hash{$id}[5]) && ($junctionID[12] >= $hash{$id}[5])) || (($junctionID[11] >= $hash{$id}[4]) && ($junctionID[12] <= $hash{$id}[5]))) {
                                  if($hash{$id}[0] ~~ @name) { }
		                  else {
			              print OUT $junctionID[0]; print OUT "\t"; print OUT $junctionID[1]; print OUT "\t"; print OUT $junctionID[2]; print OUT "\t";
                                      print OUT $hash{$id}[0]; print OUT "\n";
			              push @name, $hash{$id}[0];
                                      $tagg=1;
                                  }
                              }
                        }
       }
   if($tagg==0) {
                   print OUT $junctionID[0]; print OUT "\t"; print OUT $junctionID[1]; print OUT "\t"; print OUT $junctionID[2]; print OUT "\t";
                   print OUT "NONE"; print OUT "\n";
   }
}    
}
              
                  
                
