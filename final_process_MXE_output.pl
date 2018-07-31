#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
#data_file_1=MXE_final; $data_file2=out1; $data_file_3=out2; $data_file_4=$condition_num; $data_file_5=$ctrl_num;  

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:final_process_MXE_output.pl <data_file_1> <data_file_2> <data_file_3> <data_file_4> <data_file_5> \n";

}

my ($data_file_1, $data_file_2, $data_file_3, $condition_num, $ctrl_num)  = @ARGV;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") || die "can't open OUT1 file";
open(OUT2, ">$data_file_3") || die "can't open OUT2 file";

my %hash;
my %readnum;
my %ratio;
my $as;
my $chr;
my $strand;
my $start;
my $end;
my @id=();
my @unique_id=();
my @sorted_id=();
my @pvalue=();
my @qvalue=();
my @array;
my $pmin;
my $qmin;
my $indicator=0;
my $lineid;
my $temp;
my $count;
my $i;
my $j;
my $s;
my $l;
my $m;
my $sum;
my $sum1;
my $sum2;
my $gene;

while(<IN1>) {
      chomp;
      if($_ !~ /raw_count/) {
              if($indicator == 1) {
		      @array=split(/\s+/,$_);
		      $gene = $array[0];
		      $as = $array[1];
		      $lineid = 1;
		      $temp = $array[2];
		      $count = ($temp =~ s/_/_/g);
		      $i = $count/2;
		      $chr = $array[11]; 
		      $strand = $array[15];
                      @{ $hash{$lineid} } = @array;
		      push @id, $array[12];
		      push @id, $array[13];
		      push @pvalue, $array[6];
		      push @qvalue, $array[7];
		      $lineid=$lineid+1;
                      $indicator = $indicator + 1;
	      }
	      else {
		      @array=split(/\s+/,$_);
		      if ($array[1] eq $as) {
			      @{ $hash{$lineid} } = @array;
			      $lineid = $lineid + 1;
			      push @id, $array[12];
			      push @id, $array[13];
			      push @pvalue, $array[6];
			      push @qvalue, $array[7];
		      }
		      else { 
				      %readnum = ();
				      %ratio = ();
				      for($j=1; $j<=$i; $j++) {
				           for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
                                              push @{$readnum{$j}}, ${$hash{$j}}[16+$s]+${$hash{$j+$i}}[16+$s];   
					   }
					   # print "@{$readnum{$j}}\n";
		                      }
				       for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
					       $sum = 0;
                                                 for($j=1; $j<=$i; $j++) {
							 $sum = $sum + $readnum{$j}[$s];
						 }
						 for($j=1; $j<=$i; $j++) {
							 if($sum != 0) {
                                                           push @{$ratio{$j}}, $readnum{$j}[$s]/$sum;
					                 }
							 else {
					                   push @{$ratio{$j}}, "INF";
						        }
				               }
				       }
				       #for($j=1; $j<=$i; $j++) {
					       #print "@{$ratio{$j}}\n";
					       #}
				      @unique_id = uniq @id;
                                      @sorted_id = sort { $a <=> $b } @unique_id;
				      #print "@sorted_id\n";
				      $pmin = min @pvalue;
				      s/NA/1/g for @qvalue;
				      $qmin = min @qvalue;
                                      print OUT2 $as; print OUT2 "\t"; 
				      print OUT2 $chr; print OUT2 "_"; print OUT2 $strand; print OUT2 "_"; print OUT2 $sorted_id[0]; print OUT2 "_";
                                           for($l=1; $l<=$#sorted_id-2; $l=$l+2) {
                                                   print OUT2 $sorted_id[$l]+1; print OUT2 "-"; print OUT2 $sorted_id[$l+1]; print OUT2 "_";
                                           }
                                      print OUT2 $sorted_id[-1]; print OUT2 "\n";

				      print OUT1 $gene; print OUT1 "\t"; print OUT1 $chr; print OUT1 "_"; print OUT1 $strand; print OUT1 "_"; print OUT1 $sorted_id[0]; print OUT1 "_"; 
				           for($l=1; $l<=$#sorted_id-2; $l=$l+2) {
						   print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 "_";
					   }
				      print OUT1 $sorted_id[-1]; print OUT1 "\t"; 
			              print OUT1 $chr; print OUT1 "\t"; print OUT1 $strand; print OUT1 "\t"; print OUT1 $sorted_id[0]; print OUT1 "\t";
				           for($l=1; $l<$#sorted_id-2; $l=$l+2) {
						   print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 ";";
					   }
					   print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 "\t";
			              print OUT1 $sorted_id[-1]; print OUT1 "\t"; 
				      print OUT1 $pmin; print OUT1 "\t"; print OUT1 $qmin; print OUT1 "\t";
				      #print OUT1 join(';', @pvalue); print OUT1 "\t"; print OUT1 join(';', @qvalue); print OUT1 "\t"; 
				      for($j=1; $j<=$i-1; $j++) {
					      $sum1 = 0;
					      $sum2 = 0;
					      if ( grep( /^INF$/, @{$ratio{$j}} ) ) { print OUT1 "INF"; print OUT1 "\n";}
					      else {
			                             for($m=0; $m<$condition_num; $m++) {
			                                  $sum1 = $sum1 + $ratio{$j}[$m];
				                     }
					             for($m=0; $m<$ctrl_num; $m++) {
					                  $sum2 = $sum2 + $ratio{$j}[$condition_num+$m];
		                                     }			      
                                                    print OUT1 $sum1/$condition_num - $sum2/$ctrl_num; print OUT1 ";"
				                   }
					   }
				              $sum1 = 0;
                                              $sum2 = 0;
					      if ( grep( /^INF$/, @{$ratio{$j}} ) ) { print OUT1 "INF"; print OUT1 "\n";}
					      else {
                                                     for($m=0; $m<$condition_num; $m++) {
                                                          $sum1 = $sum1 + $ratio{$j}[$m];
						     }
						     for($m=0; $m<$ctrl_num; $m++) {
                                                          $sum2 = $sum2 + $ratio{$j}[$condition_num+$m];
                                                     }
                                                   print OUT1 $sum1/$condition_num - $sum2/$ctrl_num; print OUT1 "\n";
					      }
			              $as = $array[1];
				      $gene = $array[0];
			              %hash = ();
			              $lineid = 1;
			              $temp = $array[2];
			              $count = ($temp =~ s/_/_/g);
			              $i = $count/2;
				      $chr = $array[11];
                                      $strand = $array[15];
			              @{ $hash{$lineid} } = @array;
			              $lineid=$lineid+1;
				      @id = ();
				      @pvalue = ();
				      @qvalue = ();
				      push @id, $array[12];
				      push @id, $array[13];
				      push @pvalue, $array[6];
				      push @qvalue, $array[7];
		      }
		      $indicator = $indicator + 1;
	      }
      }
      else {
	      @array=split(/\s+/,$_);
	      print OUT1 "Gene"; print OUT1 "\t"; print OUT1 "AS_event_ID"; print OUT1 "\t"; print OUT1 "chromosome"; print OUT1 "\t"; print OUT1 "strand"; print OUT1 "\t"; print OUT1 "upstream_exon_end_coor"; print OUT1 "\t"; print OUT1 "MXE_exon_coordinates"; print OUT1 "\t"; print OUT1 "downstream_exon_start_coor"; print OUT1 "\t" ; print OUT1 "pvalue"; print OUT1 "\t"; print OUT1 "qvalue"; print OUT1 "\t"; print OUT1 "deltaPSI_"; print OUT1 $array[8]; print OUT1 "-"; print OUT1 $array[9]; print OUT1 "\n";
	      $indicator = $indicator+1;}
      #print $indicator; print "\n";
}

%readnum = ();
%ratio = ();
for($j=1; $j<=$i; $j++) {
     for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
             push @{$readnum{$j}}, ${$hash{$j}}[16+$s]+${$hash{$j+$i}}[16+$s];
     }
     #print "@{$readnum{$j}}\n";
}
for($s=0; $s<=$ctrl_num+$condition_num-1; $s++) {
                    $sum = 0;
                    for($j=1; $j<=$i; $j++) {
                               $sum = $sum + $readnum{$j}[$s];
                    }
                    for($j=1; $j<=$i; $j++) {
			    if($sum != 0) {
                               push @{$ratio{$j}}, $readnum{$j}[$s]/$sum;
                            }
			    else {
				    push @{$ratio{$j}}, "INF";
			    }
		    }
}
#for($j=1; $j<=$i; $j++) {
	#print "@{$ratio{$j}}\n";
	#}
@unique_id = uniq @id;
@sorted_id = sort { $a <=> $b } @unique_id;
#print "@id\n";
#print "@sorted_id\n";
$pmin = min @pvalue;
s/NA/1/g for @qvalue;
$qmin = min @qvalue;
print OUT2 $as; print OUT2 "\t";
print OUT2 $chr; print OUT2 "_"; print OUT2 $strand; print OUT2 "_"; print OUT2 $sorted_id[0]; print OUT2 "_";
     for($l=1; $l<=$#sorted_id-2; $l=$l+2) {
             print OUT2 $sorted_id[$l]+1; print OUT2 "-"; print OUT2 $sorted_id[$l+1]; print OUT2 "_";
     }
print OUT2 $sorted_id[-1]; print OUT2 "\n";

print OUT1 $gene; print OUT1 "\t"; print OUT1 $chr; print OUT1 "_"; print OUT1 $strand; print OUT1 "_"; print OUT1 $sorted_id[0]; print OUT1 "_"; 
     for($l=1; $l<=$#sorted_id-2; $l=$l+2) {
             print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 "_";
     }
print OUT1 $sorted_id[-1]; print OUT1 "\t"; 
print OUT1 $chr; print OUT1 "\t"; print OUT1 $strand; print OUT1 "\t"; print OUT1 $sorted_id[0]; print OUT1 "\t";
     for($l=1; $l<$#sorted_id-2; $l=$l+2) {
             print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 ";";
     }
     print OUT1 $sorted_id[$l]+1; print OUT1 "-"; print OUT1 $sorted_id[$l+1]; print OUT1 "\t";
print OUT1 $sorted_id[-1]; print OUT1 "\t";
print OUT1 $pmin; print OUT1 "\t"; print OUT1 $qmin; print OUT1 "\t"; 
#print OUT1 join(';', @pvalue); print OUT1 "\t"; print OUT1 join(';', @qvalue); print OUT1 "\t";
for($j=1; $j<=$i-1; $j++) {
   $sum1 = 0;
   $sum2 = 0;
   if ( grep( /^INF$/, @{$ratio{$j}} ) ) { print OUT1 "INF"; print OUT1 "\n";}
   else {
	   for($m=0; $m<$condition_num; $m++) {
                       $sum1 = $sum1 + $ratio{$j}[$m];
           }
           for($m=0; $m<$ctrl_num; $m++) {
                       $sum2 = $sum2 + $ratio{$j}[$condition_num+$m];
          }
           print OUT1 $sum1/$condition_num - $sum2/$ctrl_num; print OUT1 ";"
        }
 }
$sum1 = 0;
$sum2 = 0;
if ( grep( /^INF$/, @{$ratio{$j}} ) ) { print OUT1 "INF"; print OUT1 "\n";}
else {
for($m=0; $m<$condition_num; $m++) {
$sum1 = $sum1 + $ratio{$j}[$m];
}
for($m=0; $m<$ctrl_num; $m++) {
$sum2 = $sum2 + $ratio{$j}[$condition_num+$m];
}
print OUT1 $sum1/$condition_num - $sum2/$ctrl_num; print OUT1 "\n";
}



