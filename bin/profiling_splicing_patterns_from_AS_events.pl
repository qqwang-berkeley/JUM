#!/usr/bin/perl
use strict;
use warnings;
#data_file_1=temp_long_intron_retention_junction_coordinate_with_read_num.txt; date_file_2=long_intron_retention_with_linear_fitting.txt; 

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:profiling_splicing_patterns_from_AS_events.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

my %junction; #$junction{junction_id}=[# of existence in all AS events];
#my %pattern; #$pattern{number}=[AS_event1, AS_event2, ..., ]
#my %metric;  #$metric{number}=[#_of_existence_AS_1, #_of_existence_AS_2, ... , ]

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") or die "can't open input2 file: $!";
open(OUT2, ">$data_file_3") or die "can't open input3 file: $!";

while(<IN1>) {
      chomp;
      my @array=split(/\s+/,$_);
       if(exists $junction{$array[2]}) {
          $junction{$array[2]}[0]=$junction{$array[2]}[0]+1;
       }
       else {
          $junction{$array[2]}[0]=1;
       }
     my $temp=$array[0];
     my @name=split(/:/,$array[0]);
     print OUT1 $name[0]; print OUT1 "\t"; print OUT1 $temp; print OUT1 "\t"; print OUT1 $array[1]; print OUT1 "\t"; print OUT1 $array[2]; print OUT1 "\n"; 
}

close IN1; 

foreach my $i (keys %junction) {
       print OUT2 $i; print OUT2 "\t"; print OUT2 $junction{$i}[0]; print OUT2 "\n";
}
 #      my $outlier_indicator=0;
       #my $temp="none";
       #my @x=(1 .. $junction{$i}[0]);
       #print @x; print "\n";
  #     my @y = @{ $junction{$i} }[1 .. $junction{$i}[0]];
       #print @y; print "\n";
       #$lineFit = Statistics::LineFit->new();
       #$lineFit->setData (\@x, \@y) or die "Invalid data";
       #my $rSquared = $lineFit->rSquared();
       #($tStatIntercept, $tStatSlope) = $lineFit->tStatistics();
       #my $prob_t_intercept=Statistics::Distributions::tprob($junction{$i}[0]-2,$tStatIntercept);
       #my $prob_t_slope=Statistics::Distributions::tprob($junction{$i}[0]-2,$tStatSlope);
       #my $pval = dagostino_k_square_test (\@y);
   #    my $stat = Statistics::Descriptive::Full->new();
    #   $stat->add_data(\@y);
    #   my $q1 = $stat->quantile(1);
    #   my $q3 = $stat->quantile(3);
       #my $iqr = $q3-$q1;
   #    my $q0 = $stat->quantile(0);
   #    my $q4 = $stat->quantile(4);
   #    my $median = $stat->median();
       #for($j=1; $j<=$junction{$i}[0]; $j++) {
        #    if(($junction{$i}[$j] < $q1-2*$iqr) || ($junction{$i}[$j] > $q3+2*$iqr)) {
         #      $outlier_indicator=1;  $temp=$junction{$i}[$j];
          #     last;
           # }
       #}
    #   if (($q0 < $median/4) || ($q1 < $q3/3) || ($q4 > $median*4)) {
    #         $outlier_indicator=1;
    #   }              
    #   print OUT $i; print OUT "\t"; print OUT $outlier_indicator; print OUT "\t"; # print OUT $q0; print OUT "\t"; print OUT $median; print OUT "\t"; print OUT $q1; print OUT "\t"; print OUT $q3; print OUT "\t";
    #      if ( grep( /^0$/, @{ $junction{$i} } ) ) {       #this is to determine if there is any zero readout for read mapped to the intron. If yes, not valid intron retention and will be marked as 1, otherwise 0
    #        print OUT "1"; print OUT "\n";
    #      }
    #      else { print OUT "0"; print OUT "\n";
    #      }
#}
                  
       
