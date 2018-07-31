#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
#data_file_1=MXE_final; $data_file2=out1; $data_file_3=out2;  

my $len=scalar(@ARGV);
if ($len < 1) {
    die "Usage:final_process_A5SS_output.pl <data_file_1> <data_file_2> <data_file_3>\n";

}

my ($data_file_1, $data_file_2, $data_file_3)  = @ARGV;

open(IN1, "<$data_file_1") or die "can't open input1 file: $!";
open(OUT1, ">$data_file_2") || die "can't open OUT1 file";
open(OUT2, ">$data_file_3") || die "can't open OUT2 file";

my $as;
my $chr;
my $strand;
my @id=();
my @unique_id=();
my @sorted_id=();
my @pvalue=();
my @qvalue=();
my @dpsi=();
my @array;
my $pmin;
my $qmin;
my $indicator=0;
my $gene;

while(<IN1>) {
      chomp;
      if($_ !~ /raw_count/) {
              if($indicator == 1) {
		      @array=split(/\s+/,$_);
		      $gene = $array[0];
		      $as = $array[1];
		      $chr = $array[11]; 
		      $strand = $array[15];
		      push @id, $array[12];
		      push @id, $array[13];
		      push @pvalue, $array[6];
		      push @qvalue, $array[7];
		      push @dpsi, $array[-1];
                      $indicator = $indicator + 1;
	      }
	      else {
		      @array=split(/\s+/,$_);
		      if ($array[1] eq $as) {
			      push @id, $array[12];
			      push @id, $array[13];
			      push @pvalue, $array[6];
			      push @qvalue, $array[7];
			      push @dpsi, $array[-1];
		      }
		      else { 
				      @unique_id = uniq @id;
                                      @sorted_id = sort { $a <=> $b } @unique_id;
				      $pmin = min @pvalue;
				      s/NA/1/g for @qvalue;
				      $qmin = min @qvalue;
                                      print OUT2 $as; print OUT2 "\t"; 
				      print OUT2 $chr; print OUT2 "_"; print OUT2 $strand; print OUT2 "_"; print OUT2 join('_', @sorted_id); print OUT2 "\n";

				      print OUT1 $gene; print OUT1 "\t"; print OUT1 $chr; print OUT1 "_"; print OUT1 $strand; print OUT1 "_"; print OUT1 join('_', @sorted_id); print OUT1 "\t"; print OUT1 $chr; print OUT1 "\t"; print OUT1 $strand; print OUT1 "\t"; print OUT1 join(';', @sorted_id); print OUT1 "\t";
				      print OUT1 $pmin; print OUT1 "\t"; print OUT1 $qmin; print OUT1 "\t";
				      print OUT1 join(';', @dpsi); print OUT1 "\n"; 
			              
				      $as = $array[1];
				      $gene = $array[0];
				      $chr = $array[11];
                                      $strand = $array[15];
				      @id = ();
				      @pvalue = ();
				      @qvalue = ();
				      @dpsi = ();
                                      push @dpsi, $array[-1];
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
	      print OUT1 "Gene"; print OUT1 "\t"; print OUT1 "AS_event_ID"; print OUT1 "\t"; print OUT1 "chromosome"; print OUT1 "\t"; print OUT1 "strand"; print OUT1 "\t"; print OUT1 "Composite_coordinates"; print OUT1 "\t" ; print OUT1 "pvalue"; print OUT1 "\t"; print OUT1 "qvalue"; print OUT1 "\t"; print OUT1 "deltaPSI_"; print OUT1 $array[8]; print OUT1 "-"; print OUT1 $array[9]; print OUT1 "\n";
	      $indicator = $indicator+1;
      }
}

@unique_id = uniq @id; 
@sorted_id = sort { $a <=> $b } @unique_id;
$pmin = min @pvalue;
s/NA/1/g for @qvalue;
$qmin = min @qvalue;
print OUT2 $as; print OUT2 "\t"; 
print OUT2 $chr; print OUT2 "_"; print OUT2 $strand; print OUT2 "_"; print OUT2 join('_', @sorted_id); print OUT2 "\n";
print OUT1 $gene; print OUT1 "\t"; print OUT1 $chr; print OUT1 "_"; print OUT1 $strand; print OUT1 "_"; print OUT1 join('_', @sorted_id); print OUT1 "\t"; print OUT1 $chr; print OUT1 "\t"; print OUT1 $strand; print OUT1 "\t"; print OUT1 join(';', @sorted_id); print OUT1 "\t";
print OUT1 $pmin; print OUT1 "\t"; print OUT1 $qmin; print OUT1 "\t";
print OUT1 join(';', @dpsi); print OUT1 "\n";

