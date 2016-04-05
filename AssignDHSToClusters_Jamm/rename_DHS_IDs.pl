use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Embryo_stage_5/Jamm_v_1.0.6rev3/peaks/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/";


my ($original, $myID, $stage)=@ARGV;

open IN, "<$original" or print "Can't open file:$original";
open OUT, ">$myID" or print "Can't open file:$myID";
 my $line;
 

 	while ($line = <IN>)
	{
		chomp($line);

		my ($chr, $start, $end, $id, $s1, $s2, $s3, $s4, $s5, $s6) = split(/\s+/, $line);
		my $dhs_id = $chr."_".$start."_stage_".$stage;

		print OUT ("$chr\t$start\t$end\t$dhs_id\t$s1\t$s2\t$s3\t$s4\t$s5\t$s6\n");
	}
close IN;
close OUT;


