#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=output_long_intron.gff; data_file_2=XXXXXXintersect_long_intron.sam; readlength= #; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:count_reads_for_intron_retention_sam.pl <data_file_1> <data_file_2> <readlength>\n";

}

my ($data_file_1, $data_file_2, $readlength)  = @ARGV;

my %hash; #$hash{chr}=[start1, start2, start3...];

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";

#while(<IN1>) {
#	chomp;
#	my @long=split(/\s+/,$_);
#	my @longid=split(/=/,$long[8]);
#	$hash{$long[0]}{$longid[1]}[0]=$long[3];
#	$hash{$long[0]}{$longid[1]}[1]=$long[4]-$readlength;
#	$hash{$long[0]}{$longid[1]}[2]=0;
#}

#close IN1;

while(<IN2>) {
      chomp;
      my @array=split(/\s+/,$_);
      if (exists $hash{$array[0]}) {
                push @{ $hash{$array[0]} }, $array[1];
      }
      else {
	      $hash{$array[0]}[0]=$array[1];
	   }
}
close IN2;

while(<IN1>) {
	chomp;
	my @long=split(/\s+/,$_);
	my @longid=split(/=/,$long[8]);
	if(exists $hash{$long[0]}) {
		# my @sorted = sort @ { $hash{$long[0]} };
	    my @count = grep {($_ >= $long[3]) && ($_ <= $long[4]-$readlength)} @{ $hash{$long[0]} };
	    my $num = scalar @count;
	    print $longid[1]; print "\t"; print $num; print "\n";
        }
	else {
	    print $longid[1]; print "\t"; print "0"; print "\n";
        }
}	

close IN1;

