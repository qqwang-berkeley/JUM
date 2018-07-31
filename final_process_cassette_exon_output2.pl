#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=Index1_Aligned.out.spanning_junction_reads_bed12.bed; data_file_2=file2 and data_file_3=output;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage: perl final_process_cassette_exon_output2 <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %hash; #$hash{chr}{start}{end}=[leftoverhang,rightoverhang];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(OUT, ">$data_file_3") || die "can't open OUT file";
my $i;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[0]}[0]=$array[1];
}

while(<IN2>) {
	chomp;
	if($_ !~ /raw_count/) {
	   my @temp=split(/\s+/,$_);
	   $i=$hash{$temp[1]}[0];
	   $temp[1]=$i;
	   print OUT join( "\t", @temp ), "\n";
        }
	else {
		print OUT $_; print OUT "\n";
	}
}


              
                  
                
