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

my %hash; #$hash{AS_event_id}=[gene1_#, gene2_#, gene3_#, ...];
my %AS; #$AS{AS_event_id}=[gene1,gene2,gene3,...];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") || die "can't open OUT file: $!";
open(OUT2, ">$data_file_3") || die "can't open OUT file";

my @array;
my $indi;
my $l;
my @temp;
my $maxvalue;
my @maxindex;

while(<IN1>) {
      chomp;
      @array=split(/\s+/,$_);
      if(exists $AS{$array[0]}) {
	      $indi=0;
	      for($l=1;$l<=$AS{$array[0]}[0];$l++) {
		      if($AS{$array[0]}[$l] eq $array[3]) {
			      $hash{$array[0]}[$l] = $hash{$array[0]}[$l]+1;
		              $indi=1;
			      last;
		      }
              }
	      if($indi==0) {
		      $AS{$array[0]}[0]=$AS{$array[0]}[0]+1;
		      $AS{$array[0]}[$AS{$array[0]}[0]]=$array[3];
		      $hash{$array[0]}[0]=$hash{$array[0]}[0]+1;
		      $hash{$array[0]}[$hash{$array[0]}[0]]=1;
	      }
      }
      else {
	      $AS{$array[0]}[0]=1;
	      $AS{$array[0]}[1]=$array[3];
	      $hash{$array[0]}[0]=1;
	      $hash{$array[0]}[1]=1;
      }
}
close IN1;

foreach my $n (keys %AS) {
	if($AS{$n}[0] == 1) {
		print OUT1 $n; print OUT1 "\t"; print OUT1 $AS{$n}[1]; print OUT1 "\n";
	}
	else {
              @temp = @{$hash{$n}}[1 .. $hash{$n}[0]];
              $maxvalue = max @temp;  
              @maxindex = grep $temp[$_] eq $maxvalue , 0 .. $#temp;
	      if($#maxindex == 0) {
		print OUT1 $n; print OUT1 "\t"; print OUT1 $AS{$n}[$maxindex[0]+1]; print OUT1 "\n";
	      }
	      else {
	        print OUT1 $n; print OUT1 "\t"; print OUT1 $AS{$n}[$maxindex[0]+1]; print OUT1 "\n";
		for($l=1; $l<=$#maxindex; $l++) {
			print OUT2 $n; print OUT2 "\t"; print OUT2 $AS{$n}[$maxindex[$l]+1]; print OUT2 "\n";
		}
	      }
        }
}
      
                  
                
