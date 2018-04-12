#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=UNION_junc_coor_with_junction_ID_more_than_5_reads_overlap_from_all_samples_formatted_with_junction_length_with_overhang_union_fromm_all_samples.txt; data_file_2=output_long_intron.gff; data_file_3=output_short_intron.gff; readlength= #; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:extract_intron_retention_event_coordinate_from_conditional_union_splicing_junction_gff.pl <data_file_1> <data_file_2> <data_file_3> <readlength>\n";
}

my ($data_file_1, $data_file_2, $data_file_3, $readlength)  = @ARGV;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") || die "can't open OUT file";
open(OUT2, ">$data_file_3") || die "can't open OUT file";
#open(OUT3, ">$data_file_4") || die "can't open OUT file";

my $exon_intron_start;
my $exon_intron_end;
my $intron_exon_start;
my $intron_exon_end;
my $intron_retention_start;
my $intron_retention_end;
my $short_intron_left_ext;
my $short_intron_right_ext;
my $num;

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
#####################################
      $exon_intron_start=0;
      $exon_intron_end=0;
      $intron_exon_start=0;
      $intron_exon_end=0;
      $intron_retention_start=0;
      $intron_retention_end=0;
#####################################
      if($array[4] >= $readlength-2 ) {
                  $exon_intron_start=$array[2]-$array[6];
		      if($exon_intron_start == 0) {
			      $exon_intron_start = $exon_intron_start+1;
		      }
                  $exon_intron_end=$array[2]+$readlength-2;
                  $intron_exon_start=$array[3]-($readlength-2)+1;
                  $intron_exon_end=$array[3]+$array[7]+1;
      
            print OUT1 $array[0]; print OUT1 "\t"; print OUT1 "intron_retention"; print OUT1 "\t"; print OUT1 "retention"; print OUT1 "\t";print OUT1 $exon_intron_start; print OUT1 "\t"; print OUT1 $exon_intron_end; print OUT1 "\t"; print OUT1 "."; print OUT1 "\t"; if (($array[1] ne "+") && ($array[1] ne "-")) {print OUT1 ".";} else {print OUT1 $array[1];}  print OUT1 "\t"; print OUT1 "."; print OUT1 "\t"; print OUT1 "junc_id"; print OUT1 "="; print OUT1 $array[5]; print OUT1 "-"; print OUT1 "1"; print OUT1 "\n"; 
            print OUT1 $array[0]; print OUT1 "\t"; print OUT1 "intron_retention"; print OUT1 "\t"; print OUT1 "retention"; print OUT1 "\t";print OUT1 $intron_exon_start; print OUT1 "\t"; print OUT1 $intron_exon_end; print OUT1 "\t"; print OUT1 "."; print OUT1 "\t"; if (($array[1] ne "+") && ($array[1] ne "-")) {print OUT1 ".";} else {print OUT1 $array[1];} print OUT1 "\t"; print OUT1 "."; print OUT1 "\t"; print OUT1 "junc_id"; print OUT1 "="; print OUT1 $array[5]; print OUT1 "-"; print OUT1 "2"; print OUT1 "\n";
            
           # for($num=0; $num<=$array[4]-1; $num++) {
            #       print OUT3 $array[0]; print OUT3 "\t"; print OUT3 $array[2]+$num; print OUT3 "\t"; print OUT3 $array[2]+$num+1; print OUT3 "\t"; print OUT3 $array[1]; print OUT3 "\t"; print OUT3 $array[5]; print OUT3 "\n";
          # } 
      }

      else {   $short_intron_left_ext= int(($readlength-$array[4])/2 + 0.5);
               $short_intron_right_ext=int(($readlength-$array[4])/2 + 0.5);
                   if(($array[6] >= $short_intron_left_ext) && ($array[7] >= $short_intron_right_ext)) {
                      $intron_retention_start=$array[2]-$short_intron_left_ext;
                      $intron_retention_end=$array[3]+$short_intron_right_ext+1;
                   }
                   elsif(($array[6] < $short_intron_left_ext) && ($array[7] >= $short_intron_right_ext+$short_intron_left_ext-$array[6])) {
                      $intron_retention_start=$array[2]-$array[6];
                      $intron_retention_end=$array[3]+$short_intron_right_ext+$short_intron_left_ext-$array[6]+1;
                   }
                   elsif(($array[7] < $short_intron_right_ext) && ($array[6] >= $short_intron_left_ext+$short_intron_right_ext-$array[7])) {
                      $intron_retention_start=$array[2]-($short_intron_left_ext+$short_intron_right_ext-$array[7]);
                      $intron_retention_end=$array[3]+$array[7]+1;
                   }
                   else {
                      $intron_retention_start=$array[2]-$array[6];
                      $intron_retention_end=$array[3]+$array[7]+1;
                  }
               
            print OUT2 $array[0]; print OUT2 "\t"; print OUT2 "intron_retention"; print OUT2 "\t"; print OUT2 "retention"; print OUT2 "\t";print OUT2 $intron_retention_start; print OUT2 "\t"; print OUT2 $intron_retention_end; print OUT2 "\t"; print OUT2 "."; print OUT2 "\t"; if (($array[1] ne "+") && ($array[1] ne "-")) {print OUT2 ".";} else {print OUT2 $array[1];} print OUT2 "\t"; print OUT2 "."; print OUT2 "\t"; print OUT2 "junc_id"; print OUT2 "="; print OUT2 $array[5]; print OUT2 "\n";
     }
}
            
close IN1;




