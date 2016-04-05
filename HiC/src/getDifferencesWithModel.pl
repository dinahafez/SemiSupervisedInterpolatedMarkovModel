use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $i = $ARGV[0];
my $j=$ARGV[1];

my $overlapFile = "/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/stage1_12_ovrlap_hic.bed";
my $resultDrive = "/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/stage1_14_intersect_partnerTogene_stage1.bed.filtered.uniq.bed"; 
#model_uniq_genes.bed";


#"/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2_pos_".$i."/all_pos_".$j."_uniq_genes.bed";
my $output = "/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/stage14_try";

#"/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2_pos_".$i."/all_pos_".$j."_overlap_HiC.bed";
my %overlap;	
#read overlap
open IN, "<$overlapFile" or die "Can not open file :$overlapFile";
 
while (my $line = <IN> )
{
	
	chomp($line);
	
	#print ("$line\n");
	my @fields = split(/\s+/,$line);
	push(@{$overlap{$fields[3]}}, $fields[1]-1);
}
close IN;

open IN, "<$resultDrive" or die "Can not open file :$resultDrive";
open OUT, ">$output" or die "Can not open file :$output";
 
while (my $line = <IN> )
{
		
	chomp($line);
	
	my @fields = split(/\s+/,$line);
	

	if (grep(/^$fields[1]/i, @{$overlap{$fields[3]}}))
	{
		print OUT ("$line\n");
	}
	else
	{
		#print OUT ("$line\n");
	}
}
close IN;
close OUT;

	
 		
