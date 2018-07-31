#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=UNION_junc_coor_with_junction_ID_more_than_5_reads_overlap_from_all_samples_formatted_with_junction_length.txt; data_file_2=Index1_2nd_pass_junction_counts_more_than_five_in_all_samples_with_both_overhangs.txt; data_file_3=Index2_2nd_pass_junction_counts_more_than_five_in_all_samples_with_both_overhangs.txt; data_file_4=Index3_2nd_pass_junction_counts_more_than_five_in_all_samples_with_both_overhangs.txt; data_file_4=...; data_file_5=...; data_file_n=OUTPUT 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:extract_intron_retention_event_coordinate_from_conditional_union_samples.pl <input_1> <input_2> ... <output_file>\n";
}

my $mark;
my $loop;
my $o=1;
my %hash; #$hash{junction_ID}=[$chr,$strand, $start, $end, $junction_span,$left_overhang, $right_overhang];
my @file;
$file[0]=0;
my $filehandle;

for($mark=1;$mark<=$len-1;$mark++) {
        open($file[$mark], "<$ARGV[$mark-1]") || die "can't open intput file $mark";
}

open(OUT, ">$ARGV[$len-1]") || die "can't open OUT file";

$filehandle=$file[1];

while(<$filehandle>) {
      chomp;
      my @array=split(/\s+/,$_);
      $hash{$array[5]}[0]=$array[0];
      $hash{$array[5]}[1]=$array[1];
      $hash{$array[5]}[2]=$array[2];
      $hash{$array[5]}[3]=$array[3];
      $hash{$array[5]}[4]=$array[4];
      $hash{$array[5]}[5]=0;
      $hash{$array[5]}[6]=0;
}

close $filehandle;

       for ($loop=2; $loop<=$len-1; $loop++) {
             $filehandle = $file[$loop];
             while(<$filehandle>) {
                     chomp;
                     my @Xtemp=split(/\s+/,$_);
                     if($Xtemp[8] > $hash{$Xtemp[4]}[5]) {
                             $hash{$Xtemp[4]}[5]=$Xtemp[8];
                     }
                     if($Xtemp[9] > $hash{$Xtemp[4]}[6]) {
                             $hash{$Xtemp[4]}[6]=$Xtemp[9];
                     }

            }
             close $filehandle;
      }


foreach my $n (keys %hash) {
          print OUT $hash{$n}[0]; print OUT "\t"; print OUT $hash{$n}[1]; print OUT "\t"; print OUT $hash{$n}[2]; print OUT "\t"; print OUT $hash{$n}[3]; print OUT "\t"; print OUT $hash{$n}[4]; print OUT "\t"; print OUT $n; print OUT "\t"; print OUT $hash{$n}[5]; print OUT "\t"; print OUT $hash{$n}[6];print OUT "\n";
}
