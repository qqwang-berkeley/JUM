#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);
#data_file_1=non_valid_cassette_exon_list.txt; date_file_2=Valid_total_MXE_list.txt; data_file_3=total_A5SS_event.txt; data_file_4=total_A3SS_event.txt; data_file_5=Valid_total_A5SS_event.txt; data_file_6=Valid_total_A3SS_event.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:screening_A5SS_A3SS_events.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4> <data_file_5> <data_file_6>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4, $data_file_5, $data_file_6)  = @ARGV;

my %A5SS; #$A5SS{ID}=[];
my %A3SS; #$A3SS{ID}=[];
open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(IN3, "<$data_file_3") or die "can't open input3 file: $!";
open(IN4, "<$data_file_4") or die "can't open input4 file: $!";
open(OUT1, ">$data_file_5") or die "can't open input5 file: $!";
open(OUT2, ">$data_file_6") or die "can't open input6 file: $!";

while(<IN3>) {
        chomp;
        $A5SS{$_}[0]=0;
}
close IN3;

while(<IN4>) {
	chomp;
        $A3SS{$_}[0]=0;
}
close IN4;

my $i;

while(<IN2>) {
	chomp;
	my @array=split(/\*/,$_);
        for($i=0; $i<=$#array; $i++) {
		if(exists $A5SS{$array[$i]}) {
			delete($A5SS{$array[$i]});
		}
		if(exists $A3SS{$array[$i]}) {
			delete($A3SS{$array[$i]});
		}
	}
}
close IN2;

while(<IN1>) {
	chomp;
	my @line=split(/\s+/,$_);
	my @reline=split(/\*/,$line[4]);
        if(($line[3] eq "+") || ($line[3] eq "0")) {
	     $A5SS{$reline[0]}[0]=0;
	     $A3SS{$reline[1]}[0]=0;
        }
        if($line[3] eq "-") {
	     $A3SS{$reline[0]}[0]=0;
	     $A5SS{$reline[1]}[0]=0;
        }
}

close IN1; 

foreach my $l (keys %A5SS) {
        print OUT1 $l; print OUT1 "\n";
}

foreach my $m (keys %A3SS) {
	print OUT2 $m; print OUT2 "\n";
}

       
