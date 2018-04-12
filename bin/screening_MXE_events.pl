#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);
#data_file_1=total_A5SS_event.txt; date_file_2=total_A3SS_event.txt; data_file_3=reconstruct_splicing_pattern_input_1.txt ; data_file_3=MXE_bed; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:screening_cassette_exon_events_1.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4>\n";

}

my ($data_file_1, $data_file_2, $data_file_3, $data_file_4)  = @ARGV;

my %AS5; #$AS5{chr}{strand}{#}{ID}=[common;1st_coor;2nd_coor;3rd_coor;xxx];
my %AS3; #$AS3{chr}{strand}{#}{ID}=[common;1st_coor;2nd_coor;3rd_coor;xxx];
my %junction; #$junction{ID}=[chr,strand,start,end]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(IN2, "<$data_file_2") or die "can't open input2 file: $!";
open(IN3, "<$data_file_3") or die "can't open input3 file: $!";
open(OUT, ">$data_file_4") or die "can't open input4 file: $!";

while(<IN3>) {
       chomp;
       my @array=split(/\s+/,$_);
       $junction{$array[6]}[0]=$array[2];
       $junction{$array[6]}[1]=$array[3];
       $junction{$array[6]}[2]=$array[4];
       $junction{$array[6]}[3]=$array[5];
}
close IN3;

my $i;
my $temp;
my $first;

while(<IN1>) {
      chomp;
      my $ori=$_;
      my @line=split(/_Junction_/,$ori);
      if($line[0] == 5) {
         $first="Junction_" . $line[1];
         $AS5{$junction{$first}[0]}{$junction{$first}[1]}{$#line}{$ori}[0]=$junction{$first}[2];
             for($i=1;$i<=$#line;$i++) {
                $temp="Junction_" . $line[$i];
	        $AS5{$junction{$first}[0]}{$junction{$first}[1]}{$#line}{$ori}[$i]=$junction{$temp}[3];
	     }
      }
      if($line[0] == 3) {
         $first="Junction_" . $line[1];
         $AS3{$junction{$first}[0]}{$junction{$first}[1]}{$#line}{$ori}[0]=$junction{$first}[3];
             for($i=1;$i<=$#line;$i++) {
                $temp="Junction_" . $line[$i];
	        $AS3{$junction{$first}[0]}{$junction{$first}[1]}{$#line}{$ori}[$i]=$junction{$temp}[2];
	    }
      }
}
close IN1;

while(<IN2>) {
      chomp;
      my $ori=$_;
      my @line2=split(/_Junction_/,$ori);
      if($line2[0] == 5) {
           $first="Junction_" . $line2[1];
           $AS5{$junction{$first}[0]}{$junction{$first}[1]}{$#line2}{$ori}[0]=$junction{$first}[2];
              for($i=1;$i<=$#line2;$i++) {
                 $temp="Junction_" . $line2[$i];
	         $AS5{$junction{$first}[0]}{$junction{$first}[1]}{$#line2}{$ori}[$i]=$junction{$temp}[3];
              }
      }
      if($line2[0] == 3) {
           $first="Junction_" . $line2[1];
           $AS3{$junction{$first}[0]}{$junction{$first}[1]}{$#line2}{$ori}[0]=$junction{$first}[3];
               for($i=1;$i<=$#line2;$i++) {
                 $temp="Junction_" . $line2[$i];
	         $AS3{$junction{$first}[0]}{$junction{$first}[1]}{$#line2}{$ori}[$i]=$junction{$temp}[2];
               }
      }
}
close IN2;

my $j;
my $indicator;
my $s;
foreach my $k (keys %AS5) {
	     foreach my $l (keys %{$AS5{$k}}) {
		     foreach my $m (keys %{$AS5{$k}{$l}}) {
			     foreach my $n (keys %{$AS5{$k}{$l}{$m}}) {
				     foreach my $s (keys %{$AS3{$k}{$l}{$m}}) {
				      $indicator=0;
                                            for($j=1; $j<=$#{ $AS5{$k}{$l}{$m}{$n} }; $j++) {
						  if($j<$#{ $AS5{$k}{$l}{$m}{$n} }) {
							  # if(($AS5{$k}{$l}{$m}{$n}[$j] <= $AS3{$k}{$l}{$m}{$s}[$j]) && ($AS3{$k}{$l}{$m}{$s}[$j] <= $AS5{$k}{$l}{$m}{$n}[$j+1])) {
							   #if(($AS5{$k}{$l}{$m}{$n}[$j] < $AS3{$k}{$l}{$m}{$s}[$j]) && ($AS3{$k}{$l}{$m}{$s}[$j] < $AS5{$k}{$l}{$m}{$n}[$j+1])) {
                                                            if(($AS5{$k}{$l}{$m}{$n}[$j] < $AS3{$k}{$l}{$m}{$s}[$j]-1) && ($AS3{$k}{$l}{$m}{$s}[$j] < $AS5{$k}{$l}{$m}{$n}[$j+1]-1)) {

							    	 $indicator=$indicator+1;
						         }
						         else {
							   $indicator=-2;
							    last;
						         }
					           }
						  if($j == $#{ $AS5{$k}{$l}{$m}{$n} }) {
							  #if($AS5{$k}{$l}{$m}{$n}[$j] <= $AS3{$k}{$l}{$m}{$s}[$j]) {
							  #if($AS5{$k}{$l}{$m}{$n}[$j] < $AS3{$k}{$l}{$m}{$s}[$j]) {
                                                           if($AS5{$k}{$l}{$m}{$n}[$j] < $AS3{$k}{$l}{$m}{$s}[$j]-1) {

							    	 $indicator=$indicator+1;
							  }
							  else {
						           $indicator=-2;
						            last;
							  }

					           }
                                             }
				       if($indicator==$m) {
                                           for($j=1; $j<=$m; $j++) {
						   print OUT $k; print OUT "\t"; print OUT $AS5{$k}{$l}{$m}{$n}[$j]+1; print OUT "\t"; print OUT $AS3{$k}{$l}{$m}{$s}[$j]; print OUT "\t"; print OUT $l; print OUT "\t"; print OUT $n; print OUT "*"; print OUT $s; print OUT "\n";
					   }
				        }
				}
			}
		}
	}
}




					    



#print OUT $junction{$cassette[0]}[0]; print OUT "\t"; print OUT $junction{$cassette[1]}[3]+1; print OUT "\t"; print OUT $junction{$cassette[0]}[2]; print OUT "\t"; print OUT $junction{$cassette[0]}[1]; print OUT "\t"; print OUT $_; print OUT "\n"; 
#      }
#}

#close IN1; 

                  
       
