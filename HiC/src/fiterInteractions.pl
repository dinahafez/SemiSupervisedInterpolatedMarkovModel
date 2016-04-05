use strict;
#use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);
#use Statistics::Zscore;

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
#my $dataDrive =	"/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/data/";
#my $resultDrive = "/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2_pos_14/";
#my $resultDrive = "/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/";
#/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2/MergedDHS/5_mer/";
	my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/";


my %genes;
my %genesHiC;
readGenesInStudy();
filterInteraction();
#Count();
#To measure the overlap
#readGenesInHiC();
#filterInteraction_otherside();

#getDistance();

#read genes in the analysis
sub readGenesInStudy
{
	
	for (my $i=15; $i <= 38; $i++)
	{
		if (($i!= 19)&&($i!= 21)&&($i!= 22)&&($i!= 30)&&($i!= 31))
		{
			my $file = $dataDrive."geneClusters/".$i.".gene.cluster.uniq";
			open IN, "<$file" or die "Can't open file:$file";
			while (my $line = <IN>)
			{
				chomp ($line);
				my ($genename, $fbgn) = split(/\s+/,$line);
				$genes{$fbgn} = 1;
			}
			close IN;
			$file = $dataDrive."geneClusters/$i.gene.cluster.non.uniq";
			open IN, "<$file" or die "Can't open file:$file";
			while (my $line = <IN>)
			{
				chomp ($line);
				my ($genename, $fbgn) = split(/\s+/,$line);
				$genes{$fbgn} = 1;
			}
			close IN;
		}
	}


}

#filter Interaction file to only include those interactions for genes in the study
sub filterInteraction
{
	my $file = $dataDrive."HiC/partnerTogene_genes_intersect_DHS.bed";
	open IN, "<$file" or die "Can't open file:$file";

	my $fileout = $dataDrive."HiC/partnerTogene_genes_in_study.bed";
	open OUT, ">$fileout" or die "Can't open file:$fileout";
	while (my $line = <IN>)
	{
		
		chomp ($line);
		my ($chr, $start, $end, $genes1, $conf, $readCount) = split(/\s+/,$line);
		my @allGenes = split(/\,/,$genes1);
		
		for(my $k=0; $k<@allGenes; $k++)
		{
			my @someGenes = split(/\=/, $allGenes[$k]);
			my @genes2 = split(/\;/, $someGenes[0]);	
			
			for (my $l=0; $l < @genes2; $l++)
			{
				if (exists ($genes{$genes2[$l]}))
				{
					print OUT ("$chr\t$start\t$end\t$genes2[$l]\n");
				}
			}

		}
	}
	close IN;
	close OUT;
}

sub getDistance
{
	my $file = $dataDrive."sextonCell_95_partnerTogene.bed";
	open IN, "<$file" or die "Can't open file:$file";

	my $fileout = $resultDrive."distance.bed";
	open OUT, ">$fileout" or die "Can't open file:$fileout";
	my $count =0;
	while (my $line = <IN>)
	{
		
		chomp ($line);
		my ($chr, $start, $end, $genes1, $conf, $readCount) = split(/\s+/,$line);
		my @allGenes = split(/\,/,$genes1);
		
		for(my $k=0; $k<@allGenes; $k++)
		{
			my @someGenes = split(/\=/, $allGenes[$k]);
			my @genes2 = split(/\;/, $someGenes[0]);	
			my ($chrg, $coord) = split(/\:/,$someGenes[1]);
			my ($startg, $endg ) = split (/\-/,$coord);
			my $mid = (($end - $start)/2) + $start;			
			my $distance = abs($endg - $mid);
			if ($chrg eq $chr)
			{			
				print OUT ("$distance\n");
			}
			else
			{$count ++;}

		}
	}
	close IN;
	close OUT;
	print ("No of interactions that exist on between different chr = $count\n");
}


sub Count
{
	my %HicGenes;
	my %HicFrag;

	my $file = $dataDrive."sextonCell_95_partnerTogene.bed";
	open IN, "<$file" or die "Can't open file:$file";

	
	while (my $line = <IN>)
	{
		
		chomp ($line);
		my ($chr, $start, $end, $genes1, $conf, $readCount) = split(/\s+/,$line);
		my @allGenes = split(/\,/,$genes1);
		
		my $id = $chr."_".$start;

		$HicFrag{$id} = 1;
		for(my $k=0; $k<@allGenes; $k++)
		{
			my @someGenes = split(/\=/, $allGenes[$k]);
			my @genes2 = split(/\;/, $someGenes[0]);	
			for (my $l=0; $l < @genes2; $l++)
			{
				
				$HicGenes{$genes2[$l]}=1;
			}
			
		}
	}
	close IN;


	my $g = keys %HicGenes;
	my $f = keys %HicFrag;

	print ("Total number of HiC Genes = $g\n");
	print ("Total number of HiC Fragments = $f\n");
}



#read genes in the analysis

sub readGenesInHiC
{
	my $file = $dataDrive."partnerTogene_genes_in_analysis.bed";#"sextonCell_95_partnerTogene.bed";
	open IN, "<$file" or die "Can't open file:$file";
	while (my $line = <IN>)
	{
		
		chomp ($line);
		my ($chr, $start, $end, $genes1, $conf, $readCount) = split(/\s+/,$line);
		#my @allGenes = split(/\,/,$genes1);
		
		#for(my $k=0; $k<@allGenes; $k++)
		#{
		#	my @someGenes = split(/\=/, $allGenes[$k]);
		#	my @genes2 = split(/\;/, $someGenes[0]);	
			
		#	for (my $l=0; $l < @genes2; $l++)
		#	{
		#		print ("$genes2[$l]\t");
		#		$genesHiC{$genes2[$l]}=1;
				
		#	}

		#}
			$genesHiC{$genes1} = 1;
	}
	close IN;
}



#filter Interaction file to only include those interactions for genes in the study
sub filterInteraction_otherside
{
	my $file = $ARGV[0]; #model_not_validated_HiC_intersect_HiC.bed ";
	open IN, "<$file" or die "Can't open file:$file";

	my $fileout =$ARGV[0].".filtered";#"model_not_validated_HiC_intersect_HiC_and_genes.bed";
	open OUT, ">$fileout" or die "Can't open file:$fileout";
	while (my $line = <IN>)
	{
		
		chomp ($line);
		my ($chr, $start, $end, $genes1) = split(/\s+/,$line);
		if (exists ($genesHiC{$genes1}))
		{
			print OUT ("$line\n");
		}
	}
	close IN;
	close OUT;

	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 ".$fileout." --uniq > ".$fileout.".uniq.bed";
	system($cmd);

}


