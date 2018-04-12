#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=Index1_2nd_pass_junction_counts_combined_intron_retention_event_list.txt; data_file_2=Index2_2nd_pass_junction_counts_combined_intron_retention_event_list.txt; data_file_3=...; data_file_4=...; data_file_5=...; data_file_6=OUTPUT; data_file_n=file_number;

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:Identify_intron_retention_event_exist_in_all_samples.pl <input_1> <input_2> ... <output_file> <file_number>\n";
}

my $mark;
my $loop;
my $o=1;
my %hash; #$hash{junction_ID}=[$index];
my @file;
$file[0]=0;
my $filehandle;
my $file_num = $ARGV[$len-1];

for($mark=1;$mark<=$len-2;$mark++) {
        open($file[$mark], "<$ARGV[$mark-1]") || die "can't open intput file $mark";
}

open(OUT, ">$ARGV[$len-2]") || die "can't open OUT file";

       for ($loop=1; $loop<=$len-2; $loop++) {
             $filehandle = $file[$loop];
             while(<$filehandle>) {
                     chomp;
                     my @Xtemp=split(/\s+/,$_);
                     if(exists $hash{$Xtemp[0]}) {
                         $hash{$Xtemp[0]}[0]=$hash{$Xtemp[0]}[0]+1;
                     }
                     else {
                         $hash{$Xtemp[0]}[0]=1;
                    }
            }
             close $filehandle;
      }


foreach my $n (keys %hash) {
             if ($hash{$n}[0] >= $file_num) {
                          print OUT $n; print OUT "\n";
            }
}
