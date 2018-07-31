#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=AS_differential.txt; data_file_2=output; input3=total file #

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:process_raw_AS_differential_output.pl <data_file_1> <data_file_2> <files>\n";

}

my ($data_file_1, $data_file_2, $files)  = @ARGV;

my %AS_structure_count; #$AS_structure_count{structure_ID}{sub_junction_ID}=[sample1_count, sample2_count, sample3_count, sample4_count, ...];
my %AS_structure_percentage; #$AS_structure_percentage{structure_ID}{sub_junction_ID}=[sample1_%, sample2_%, sample3_%, sample4_%, ...]
my %AS_structure; #$AS_structure{structure_ID}{sub_junction_ID}=[dispersion, stat, pvalue, padj, chr, start, end, strand] 
#my %pattern; #$pattern{number}=[AS_event1, AS_event2, ..., ]
#my %metric;  #$metric{number}=[#_of_existence_AS_1, #_of_existence_AS_2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT, ">$data_file_2") or die "can't open input2 file: $!";

my $i;
my $q;
my $h;
my @name;
my @percent;
my @para;

while(<IN1>) {
      chomp;
      if($_ =~ /pvalue/) {
         my @title=split(/\s+/,$_);
         my $para_len=($#title-$files-1-5)-6;
         @name = ((0) x $files);
         @percent = ((0) x $files);
         @para = ((0) x $para_len);
            for($q=$#title-$files; $q<=$#title-1; $q++) {
                $title[$q] =~ s/countData/raw_count/;
                $name[$q-($#title-$files)]=$title[$q];
                $title[$q] =~ s/raw_count/percentage_usage/;  
                $percent[$q-($#title-$files)]=$title[$q];
                #print $name[$q-($#title-$files)]; print "\t"; print $percent[$q-($#title-$files)]; print "\n";
            }
            for($h=7; $h<=$#title-$files-1-5; $h++) {
                $title[$h] =~ s/log2fold/fitting_parameter_log2fold_change/;
                $para[$h-7]=$title[$h];
                #print $para[$h-7]; print "\n";
           }         
       }
      if($_ !~ /pvalue/) {
         my @array=split(/\s+/,$_);
         $array[2] =~ s/E/J/;
	 #$AS_structure{$array[1]}{$array[2]}[0]=$array[4]; #print $AS_structure{$array[1]}{$array[2]}[0]; print "\t";
	 #$AS_structure{$array[1]}{$array[2]}[1]=$array[5]; #print $AS_structure{$array[1]}{$array[2]}[1]; print "\t";
	 #$AS_structure{$array[1]}{$array[2]}[2]=$array[6]; #print $AS_structure{$array[1]}{$array[2]}[2]; print "\t";
	 #$AS_structure{$array[1]}{$array[2]}[3]=$array[7]; #print $AS_structure{$array[1]}{$array[2]}[3]; print "\t";
	 #$AS_structure{$array[1]}{$array[2]}[4]=$array[8];
	 #$AS_structure{$array[1]}{$array[2]}[5]=$array[9];
	 #$AS_structure{$array[1]}{$array[2]}[6]=$array[10];
	 #   for($i=$#array-$files-4-5+1; $i<=$#array-4-$files; $i++) {
	 #            $AS_structure{$array[1]}{$array[2]}[$i-($#array-$files-4-5+1)+7]=$array[$i];
                     #print $AS_structure{$array[1]}{$array[2]}[$i-($#array-$files-4-5+1)+4]; print "\t"; #print $array[$i]; print "\n";
		     #   }
           #print "\n";
	 for($i=4; $i<=$#array-4-$files; $i++) {
		  $AS_structure{$array[1]}{$array[2]}[$i-4]=$array[$i];
	  }
            for($i=$#array-$files-4+1; $i<=$#array-4; $i++) {
               $AS_structure_count{$array[1]}{$array[2]}[$i-($#array-$files-4+1)]=$array[$i];
               #print $i-($#array-$files-4+1); print "\t"; print $array[$i]; print "\n"; 
               #print $i-($#array-$files-4+1)+4; print "\t"; print $array[$i-5]; print "\n";
            }
           
     }
}
close IN1;

my $j;
foreach my $id (keys %AS_structure_count) {
    my @sum = ((0) x $files);
       for($j=0; $j<=$#sum; $j++) {
             foreach my $sub (keys %{$AS_structure_count{$id}}) { 
                    $sum[$j]=$sum[$j]+$AS_structure_count{$id}{$sub}[$j];
             }
#       print $sum[$j]; print "\n";
             foreach my $jun (keys %{$AS_structure_count{$id}}) {
#                   print $AS_structure_count{$id}{$jun}[$j]; print "\t"; print $sum[$j]; print "\t"; print $AS_structure_count{$id}{$jun}[$j]/$sum[$j]; print "\n";
                 if($sum[$j] > 0) {
                    $AS_structure_percentage{$id}{$jun}[$j]=$AS_structure_count{$id}{$jun}[$j]/$sum[$j];
                    $AS_structure_percentage{$id}{$jun}[$j] = sprintf '%.2f%%', 100 * $AS_structure_percentage{$id}{$jun}[$j];
#                    print $AS_structure_percentage{$id}{$jun}[$j]; print "\n";
                 }
		 else {$AS_structure_percentage{$id}{$jun}[$j]="INF";
		       $AS_structure_percentage{$id}{$jun}[$j] = "INF";
		} 
	    }
      }
}

#foreach my $ke (keys %AS_structure_count) {
#         foreach my $su (keys %{$AS_structure_count{$ke}}) {
#                print $AS_structure_count{$ke}{$su}[0]; print "\t"; print $AS_structure_count{$ke}{$su}[1]; print "\t"; print $AS_structure_count{$ke}{$su}[2]; print "\t"; print $AS_structure_count{$ke}{$su}[3]; print "\n";
#        }
#}

#foreach my $k (keys %AS_structure_percentage) {
#         foreach my $s (keys %{$AS_structure_percentage{$k}}) {
#                print $AS_structure_percentage{$k}{$s}[0]; print "\t"; print $AS_structure_percentage{$k}{$s}[1]; print "\t"; print $AS_structure_percentage{$k}{$s}[2]; print "\t"; print $AS_structure_percentage{$k}{$s}[3]; print "\n";
#        }
#}

my $e;
my $y;
print OUT "AS_structure_ID"; print OUT "\t"; print OUT "sub_junction_ID"; print OUT "\t"; print OUT "sub-junction_dispersion_estimate"; print OUT "\t"; print OUT "LRT_statistic-full_vs_reduced"; print OUT "\t"; print OUT "LRT_p_value-full_vs_reduced"; print OUT "\t"; print OUT "BH_adjusted_p-values"; 
   for($y=0; $y<=$#para; $y++) {
       print OUT "\t"; print OUT $para[$y];
   }
print OUT "\t"; print OUT "sub_junction_chr"; print OUT "\t"; print OUT "sub_junction_start_coor"; print OUT "\t"; print OUT "sub_junction_end_coor"; print OUT "\t"; print OUT "sub_junction_size"; print OUT "\t"; print OUT "sub_junction_strand";
   for($e=0; $e<=$#name; $e++) {
       print OUT "\t"; print OUT $name[$e];
   }
   for($e=0; $e<=$#percent; $e++) {
       print OUT "\t"; print OUT $percent[$e];
};
print OUT "\n";        

my $m;
my $n;
my $p;
foreach my $k (keys %AS_structure) {
     foreach my $s (sort keys %{$AS_structure{$k}}) {
            print OUT $k; print OUT "\t"; print OUT $s;
                  for($m=0; $m <= $#{ $AS_structure{$k}{$s}}; $m++) {
                       print OUT "\t"; print OUT $AS_structure{$k}{$s}[$m]; 
                  }
                  for($n=0; $n<=$#{ $AS_structure_count{$k}{$s}}; $n++) {
                       print OUT "\t"; print OUT $AS_structure_count{$k}{$s}[$n]; 
                  }
                  for($p=0; $p<=$#{ $AS_structure_percentage{$k}{$s}}; $p++) {            
                       print OUT "\t"; print OUT $AS_structure_percentage{$k}{$s}[$p];
                 }
           print OUT "\n";
    }
}
