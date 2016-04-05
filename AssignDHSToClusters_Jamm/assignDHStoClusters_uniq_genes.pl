use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

#This file assigns DHSs to gene clusters based on several concepts. 
#Base condition: 3 different ways
#1. Each DHS is assigned to its closest gene
#2. Each gene is assinged to its closest DHS
#3. All DHSs in +/-50kb around a gene are assigned to it
#It then gets all DHS in the specified window and generates all files needed to run the McEnahncer model
# McEnhancer is run in two stages, first in considers unique genes, then it gets all common genes
# The last couple of functions generate files for DHSs controlling ubiquitous expression   

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/";
my $codeVersion = "BaseLine_assign_DHS_to_closest_gene/";

# intersectBed -a idr.2.allStages.filtered.peaks.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -u > DHS_overlap_TSS_2bp.bed
# intersectBed -a idr.2.allStages.filtered.peaks.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -v > DHS_no_overlap_TSS_2bp.bed


my %flymap;
addID();
getGeneTSS();
getTSSperCluster();
mapGeneIDs();
calculateDHSLength();
splitDHSIntoClusters_specifiedGenes();
getmidpointDHS("/data/ohler/Dina/Drosophila/Research/data/DHS_row_data//DHS_no_overlap_TSS_2bp.bed","/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed");
getClosestPerClusterForUnspecifiedGenes();
concatenate_ubiquitous();
getUnlabeled();

###################closest Gene -- base line
getmidpointDHS("/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed","/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed");
splitDHSIntoClusters_closest_gene();
concatenate_ubiquitous_base();

################baseline2 the closest gene
getTheClosestDHSperCluster();
#getUnlabeled_closestDHSperCluster();
concatenate_ubiquitous_closestDHSperCluster();
 
 ################baseline 3_all 50kb
all_DHS_50kb();
concatenate_ubiquitous_all_50kb();

############Get DHS TSS for classification
splitDHSIntoClusters_TSS();
concatenate_ubiquitous_TSS();

############Second stage MC
getTSSperCluster_second();
filterClusterTSS();

#These two functions are to get ubiquitous for all genes not only unique ones
splitDHSIntoClusters_specifiedGenes_2nd_stage();
getClosestPerClusterForUnspecifiedGenes_2nd_stage();
concatenate_ubiquitous_2nd_stage();

getUnlabeled_2nd_stage();

second stage TSS
splitDHSIntoClusters_TSS_all_genes();
concatenate_ubiquitous_TSS_all_genes();


for (my $i=11;$i<=39; $i++)
	{
		my $cmd = "wc -l ".$resultDrive."Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output_0.01/5_mer/".$i."_pos_count.tab";
		
		system($cmd);
	}

sub mapGeneIDs
{
	my $mapfile = $dataDrive."Genome/Gene_map_new_filtered_uniq.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr, $start, $end, $strand) = split(/\s+/, $line);
 		if (exists $flymap{$ensemble} && ($flymap{$ensemble} ne $flybaseID))
 		{
 			print ("something is wrong \n");
 			
 		}
 		else
 		{
 			$flymap{$ensemble} =  $flybaseID;
 		}
 	}	
 	close IN;

}

sub addID
{
	my @stages = (5,9,10,11,14);
	foreach (my $i =0; $i < @stages ; $i++ )
	{
		my $peakFile = $dataDrive."DHS_row_data/Jamm_peaks/Embryo_stage_".$stages[$i]."/peaks/filtered.peaks.narrowPeak";
		my $outputFile = $dataDrive."DHS_row_data/Jamm_peaks/Embryo_stage_".$stages[$i].".peaks";
		open IN, "<$peakFile" or die "Can not open file :$peakFile";
 		open OUT, ">$outputFile" or die "Can't open file : $outputFile";
 		
 		my $line ;
 		while ($line = <IN>)
 		{
 			chomp($line);
 			my ($chr, $start,$end, $id, @rest) = split(/\s+/, $line);
 			$id =~ s/\./\_/g;
 			my $new_id = $id."_stage_".$stages[$i];
 			print OUT ("$chr\t$start\t$end\t$new_id\t@rest\n");
 		}
		close IN;
		close OUT;
	}
}


sub getGeneTSS
{
	my $mapfile = $dataDrive."Genome/Gene_map_new_filtered_uniq.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	
	my $tssfile = $dataDrive."Genome/Gene_TSSs_2bp.bed";
	open OUT, ">$tssfile" or die "Can not open file :$tssfile";
	
	my $line ;
 	
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID,$chr, $start, $end, $strand) = split(/\s+/, $line);
 		my $tss_start;
 		my $tss_end;
 		if ($strand eq "1")
 		{
 			$tss_start = $start -1;
 			$tss_end = $start +1;
 		}
 		elsif ($strand eq "-1")
 		{
 			$tss_start = $end -1;
 			$tss_end = $end + 1;
 			
 		}
 		else
 		{ print ("error\n");}
 		
 		if ($tss_start < 0)
 		{
 			$tss_start = 0;
 		}
 		print OUT ("chr$chr\t$tss_start\t$tss_end\t$ensemble\t$flybaseID\t$strand\n");
 		
 	}	
 	close IN;
 	close OUT;

}
#sort -k1,1 -k2,2n -k3,3n -k4,4 Genome/Gene_TSSs_2bp.bed --uniq> Genome/Gene_TSSs_2bp_uniq.bed 


sub getTSSperCluster
{
	for (my $i=1; $i<=39;$i++)
 	{
 		my $cluster_file = $dataDrive . "/geneClusters/$i.gene.cluster.non.uniq";		
 		open IN, "<$cluster_file" or die "Can not open file :$cluster_file";
 		my %genes;
 		my $line;
 		while ($line = <IN> )
	 	{
	 		chomp($line);
	 		my ($cg_gene,$flybaseID) = split(/\s+/, $line);
	 		$genes{$flybaseID} = 1;
	 	}
	 	close IN;
	 	my $noOfGenes = keys %genes;
 		
		my $tssFile = $dataDrive. "Genome/Gene_TSSs_2bp.bed";
		open IN, "<$tssFile" or die "Can not open file :$tssFile";
		my $tss_geneFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		open OUT, ">$tss_geneFile" or die "Can't open file :$tss_geneFile";
		my $line;
		my $count=0;
		my @remove;
	 	while ($line = <IN> )
	 	{
	 		chomp($line);
	 		my ($chr,$start,$end,$gene,$flybaseID,$strand) = split(/\s+/, $line);
	 		if (exists $genes{$flybaseID})
	 		{
	 			print OUT ("$line\n");
	 			$count++;
	 			push (@remove, $flybaseID);
	 			
	 		}
	 		
	 	}
	 	my %removed;
		@removed{ @remove } = delete @genes{ @remove };
	
		my $leftKeys = keys %genes;
		
	 	if($leftKeys ==0)
	 	{
	 		print ("cluster $i is done nicely\n");
	 	}
	 	else
	 	{
	 		print ("cluster $i genes = $noOfGenes, found only is $count, still remaining$leftKeys\n");
	 		print Dumper \%genes;
	 	}
	 	close OUT;
 	}
}

#cat Embryo_stage_5.peaks Embryo_stage_9.peaks Embryo_stage_10.peaks Embryo_stage_11.peaks Embryo_stage_14.peaks > Embryo_all_stages.peaks #count=151842
#intersectBed -a Embryo_all_stages.peaks -b ../../Genome/Gene_TSSs_2bp.bed -u > DHS_overlap_TSS_2bp.bed #count=27995  --used to be 12405 
#(Before I was taking the average of the TSS for all transcripts, now, I am having each transcript separately)
#intersectBed -a Embryo_all_stages.peaks -b ../../Genome/Gene_TSSs_2bp.bed -v > DHS_no_overlap_TSS_2bp.bed #count=123847   --used to be 132581
# sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq DHS_no_overlap_TSS_2bp.bed > DHS_no_overlap_TSS_2bp_uniq.bed  #this shouldnot be done....coz if the DHS has the same coorinates but open in different stages ,each should be used separately 

sub calculateDHSLength
{
	my $overlapFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_10bp_uniq.bed";
	my $overlapStats = $resultDrive."DHS_statistics/DHS_no_overlap_TSS_10bp_new.bed";
	
 	open IN, "<$overlapFile" or die "Can not open file :$overlapFile";
 	open OUT, ">$overlapStats" or die "Can't open file : $overlapStats";
 	my %overlap;
 	while (my $line = <IN> )
 	{
 		chomp($line);
 		my ($chr, $start,$end, $id) = split(/\s+/, $line);
 		my $length = abs($end - $start);
 		print OUT ("$length\n");		
 	}
 	close IN;
 	close OUT;
}

sub getmidpointDHS
 {
 	#my $bedFile = "../../data/DHS/DHS_nooverlap_Redfly.bed";
 	#my $midFile = "../../data/DHS/DHS_nooverlap_Redfly_midpoints.bed";
 	my ($bedFile,$midFile) = @_;
 	open IN, "<$bedFile" or die "Can not open file :$bedFile";
 	open OUT, ">$midFile" or die "Can't open file : $midFile";
 	
 	my %dhs;
 	my $count = 0;
 	while (my $line = <IN> )
 	{
 			chomp($line);
 		my ($chr,$start,$end,$id) = split(/\s+/, $line);
 		my $mid = int(($start+$end) /2);
 		my $mid_end = $mid +1;
 		print OUT ("$chr\t$mid\t$mid_end\t$id\n");
 	}
 	close IN;
 	close OUT;		
 }
###############Initialization Red Fly
#intersectBed -a DHS_no_overlap_TSS_2bp.bed -b ../../CRM/redfly_all_specified_modified.bed -wb  | sort -k4,4 -k14,14 --uniq > DHS_overlap_redfly_specified_genes.bed

sub splitDHSIntoClusters_specifiedGenes
{
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_specified_genes.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
	
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$gene) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
		
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_specified_genes_used_in_initialization.bed";
	open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_specified_redfly.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($cg_name,$gene) = split(/\s+/,$line);
 		
 				
 				if(exists($gene_DHS{$gene}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 				{
	 					my $s = $gene_start{$gene}[$j];
	 					my $e = $gene_end{$gene}[$j] ;
	 					print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 					print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$gene}++;
	 			}
 			
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_specified_redfly_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_specified_redfly_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_specified_redfly_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 	close ALL;
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
		
 
}

sub getClosestPerClusterForUnspecifiedGenes
{
	#intersectBed -a DHS_no_overlap_TSS_10bp.bed -b ../../CRM/redfly_all_unspecified.bed -u  > DHS_overlap_redfly_unspecified_genes.bed
	#closestBed -a DHS_overlap_redfly_unspecified_genes_midpoints.bed -b ../../Genome/fly_TSS_10bp_uniq.bed -d -t all > DHS_overlap_redfly_unspecified_genes_closest.bed
	#read mapping data
	#flybase ID could be linked to more than one gene
	my $mapfile = $dataDrive."Genome/FlyBase_Genes_map.tab";

	my %flymap;
	open IN, "<$mapfile" or die "Can not open file :$mapfile";

	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $annotationSymbo, $symbol) = split(/\s+/, $line);
 	
 
 		push (@{$flymap{$ensemble}}, $flybaseID);
 	
 	}	
 	close IN;
	
	
	
	my %geneClusters;
	my $line;
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive."geneClusters/".$i.".gene.cluster.uniq";
		open IN, "<$clusterFile" or die "Can't open file:$clusterFile";
		while ($line = <IN>)
		{
			chomp($line);
			#my @flybase = $flymap{$line};
			push (@{$geneClusters{$i}}, $line);
		}
	}
	close IN;
	
	#to get the actual start and end
	my %dhs;
	my $dhsFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes.bed";
	 
	open IN, "<$dhsFile" or die "Can't open file:$dhsFile";
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	my $closestFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes_closest.bed";
	open IN, "<$closestFile" or die "Can't open file:$closestFile";
	#my %geneClusters_clone = dclone(%geneClusters);
	my %closestGenesDHS;
	my %GeneCount;
	my %DHSCount;
	my %dhs_gene;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$dhs_mid_start,$dhs_mid_end,$dhs_id,$chr_g,$gene_tss_s,$gene_tss_e,$gene,$score, $strand,$distance) = split(/\s+/, $line);
		
		for (my $i=1;$i<=39; $i++)
		{	
			#my $flybase = $flymap{$gene}[0];		
			my @output = grep {$gene eq $_} @{$geneClusters{$i}};
			if(@output >0)
			{
				$DHSCount{$i}{$dhs_id}=1;
				$GeneCount{$i}{$gene}++;
				push(@{$closestGenesDHS{$i}}, $dhs_id);
				$dhs_gene{$dhs_id}=$gene;
				#print ("$gene \t cluster $i\n");
			}
		}
		
	}
	
	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$GeneCount{$i}};
		my $s = keys %{$DHSCount{$i}};
		print ("cluster $i\t$c\t$s\n");
	}
	

	my $allUNSPEC = $resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_unspecified_genes_used_in_initialization.bed";
	open UNSPEC, ">$allUNSPEC" or die "Can't open file:$allUNSPEC";
	
	for (my $i=1;$i<=39; $i++)
	{
		my $index=1;
		my $outputFile = $resultDrive.$codeVersion."/ClustersDHS/".$i."_closest_dhs_unspecified_genes.bed";
		open OUT, ">$outputFile" or die "Can't open file:$outputFile";
		if(exists $closestGenesDHS{$i})
		{
			for(my $l=0; $l <@{$closestGenesDHS{$i}}; $l++)
			{
				my $id = $closestGenesDHS{$i}[$l];
				my $s = $dhs{$id}{'start'};
				my $e = $dhs{$id}{'end'};
				print OUT ("$dhs{$id}{'chr'}\t$s\t$e\t$flymap{$dhs_gene{$id}}[0]_".$id."_initial\n");
				print UNSPEC ("$dhs{$id}{'chr'}\t$dhs{$id}{'start'}\t$dhs{$id}{'end'}\t$id\tcluster_$i\t$flymap{$dhs_gene{$id}}[0]\n");
				$index++;
			}
		}
		close OUT;
	
		my $redFile = $resultDrive.$codeVersion."/ClustersDHS/".$i."_specified_redfly.bed";
		my $allFile = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly.bed";
		my $cmd = "cat $redFile $outputFile > $allFile";
		system($cmd);
		
		my $mergedFile = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly_uniq.bed";
		$cmd= "sort -k1,1 -k2,2n -k3,3n --uniq $allFile > $mergedFile";	
		system($cmd);
		

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly_uniq.fa";
		my $cmd = "twoBitToFa -bed=$mergedFile -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
	}
	
		close UNSPEC;
		
		my $cat = "cat ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_specified_genes_used_in_initialization.bed ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_unspecified_genes_used_in_initialization.bed > ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_all_used_in_initialization.bed";
		system($cat);
		
	for (my $i=1;$i<=39; $i++)
	{
		my $cmd = "wc -l ".$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly_uniq.bed";
		
		system($cmd);
	}
		
		
}

sub concatenate_ubiquitous
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_specified_redfly_uniq.bed " ;#  "_initialize_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

sub getUnlabeled
{
	
	my %noexists;
	#get DHS used in initialization
	my %DHS_Initial;
	for (my $j=1; $j<= 39; $j++)
	{
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/".$j."_specified_redfly_uniq.bed";#_initialize_redfly_uniq.bed ";
		open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
		
			chomp($line);
			my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
			my ($gene,$chr,$id,$stageW,$stage,$state ) = split(/\_/,$dhs_id);
			my $theID = $chr."_".$id."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
		close IN;	
	}
	my $k = keys %DHS_Initial;
	print ("number of DHS in initialization = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/DHS_unlabeled.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if(exists $DHS_Initial{$dhs_id})
		{
			#print ("DHS $dhs_id does  exist\n");
		}
		else{
			print OUT  ("$chr\t$start\t$end\t$dhs_id\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/DHS_unlabeled_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		#print $cmd;
		system($cmd);
	}

	for (my $i=1;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$genegc,$geneflybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			
			print OUTL ("$chr_dhs\t$s\t$e\t$geneflybase"."_".$dhs_id."_unlabeled\n");
			
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
	
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
	foreach my $key (keys %noexists)
	{
		print ("$key\n");
	}
	
}



#############TSS 
#intersectBed -a Embryo_all_stages.peaks -b ../../Genome/Gene_TSSs_2bp.bed -wb > DHS_overlap_TSS_2bp_with_gene_names.bed

sub splitDHSIntoClusters_TSS
{
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$cg_gene, $gene, $strand) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($gene,$flybase) =  split(/\s+/, $line);
 			
 				if(exists($gene_DHS{$flybase}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
	 				{
	 					my $s = $gene_start{$flybase}[$j];
	 					my $e = $gene_end{$flybase}[$j] ;
	 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$flybase"."_".$gene_DHS{$flybase}[$j]."_TSS\n");
#	 					print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$flybase}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$flybase}++;
	 			
 			}
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
}

sub concatenate_ubiquitous_TSS
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_TSS_uniq.bed ";  # _initialize_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}
###############BASE line
sub splitDHSIntoClusters_closest_gene
{
	my %dhs;
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	open IN, "<$dhs_file" or die "Can't open file:$dhs_file";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id,@rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}

	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_closest_gene.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$s1,$s2,$s3,$s4,$s5,$s6,$chr_r,$start_r,$end_r,$gene,$flybase,$strand,$distance) = split(/\s+/, $line);
		if(abs($distance) < 5000)  #DHS should be less than 5kb away from the gene
		{
			push(@{$gene_DHS{$flybase}},$dhs_id);
			my $start = $dhs{$dhs_id}{'start'};
			my $end = $dhs{$dhs_id}{'end'};
			push(@{$gene_start{$flybase}},$start);
			push(@{$gene_end{$flybase}},$end);
			push(@{$gene_chr{$flybase}},$chr);
		}
 	}
	close IN;
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
 	for (my $i=1; $i<=39; $i++)
 	{
		
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($gene, $flybase) = split(/\s+/, $line);
 			#my $gene = $flymap{$line};
 			
 			if(exists($gene_DHS{$flybase}))
 			{
 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
 				{
 					my $s = $gene_start{$flybase}[$j];
 					my $e = $gene_end{$flybase}[$j] ;
 					my $id = $flybase."_".$gene_DHS{$flybase}[$j];
 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$id\n");
 					$dhsCount{$i}++;
 					
 				}
 				$geneCount{$i}{$flybase}++;
 			}
 		}

		close IN;

		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.non.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($gene, $flybase) = split(/\s+/, $line);
 			#my $gene = $flymap{$line};
 			
 			if(exists($gene_DHS{$flybase}))
 			{
 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
 				{
 					my $s = $gene_start{$flybase}[$j];
 					my $e = $gene_end{$flybase}[$j] ;
 					my $id = $flybase."_".$gene_DHS{$flybase}[$j];
 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$id\n");
 					$dhsCount{$i}++;
 					
 				}
 				$geneCount{$i}{$flybase}++;
 			}
 		}


 		close OUT;
 		close IN;
 		
 		my $uniq = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos.uniq";
 		my $cmd = "uniq $clusterFile_fa > $uniq";
 		print ($cmd);
 		print ("\n");
 		system($cmd);
 		
 		my $sorted = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos.sorted";
 		my $cmd = "sort  -k1,1 -k2,2n -k3,3n -k4,4 --uniq $uniq > $sorted";
 		
 		print ($cmd);
 		print ("\n");
 		system($cmd);
 		
 		my $pos_out_file = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos.fa";
 		my $cmd = "twoBitToFa -bed=$sorted -noMask ../../data/Genome/Drosophila_R5.2bit $pos_out_file";
 		print ($cmd);
 		print ("\n");
 		system($cmd);
		
 	}
 	for (my $i=1;$i<=39; $i++)
	{
				
		my $pos_out_file = $resultDrive."/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos.fa";
		my $pos_out_file_single = $resultDrive."/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$i."_pos_single.fa";
		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $pos_out_file > $pos_out_file_single";
		print ($cmd);
 		print ("\n");
 		system($cmd);
		
	}
 	

}

sub concatenate_ubiquitous_base
{
	my $output = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$j."_pos.sorted ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = "rm ".$resultDrive."BaseLine_assign_DHS_to_closest_gene/ClustersDHS/".$j."_pos_single.fa";
		system($cmd);		
	}
}


##############################Initializa the closest gene
sub getTheClosestDHSperCluster 
{
 	
 	#read all DHSs
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";	
 	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	close IN;


	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed";
	for (my $i=1; $i<=39;$i++)
 	{
		
 		my $cluster_file = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq";	
 		my $output = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_uniq.bed";
 		my $cmd = "closestBed -a $cluster_file -b $dhs_file -d -t first >  $output";
 		print ("$cmd\n");
 		system($cmd);
		
 		
 	}
 	
 	for (my $i=1; $i<=39;$i++)
 	{
		
 		my $cluster_file = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";	
 		my $output = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_non_uniq.bed";
 		my $cmd = "closestBed -a $cluster_file -b $dhs_file -d -t first >  $output";
 		print ("$cmd\n");
 		system($cmd);

		my $all =$resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest.bed";
 		$cmd = "cat ".$resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_uniq.bed ".$resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_non_uniq.bed > ".$all;
		system($cmd);

		
 		
 	}

 	for (my $i=1;$i<=39; $i++)
	{
		
		my $output = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest.bed";
		my $bedIn = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_in.bed";
		open IN, "<$output" or die "Can't open file:$output";
		open OUT, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$flybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id, $dist) = split(/\s+/, $line);
			
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			if ($dist < 5000)
			{
				print OUT ("$chr_dhs\t$s\t$e\t".$flybase."_".$dhs_id."\n");
			}
			else
			{print ("dist > 5000 $flybase $dist\n");}
		}
		close OUT;
		
 	
 	
 		my $intialFileInUniq = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $bedIn  >$intialFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $initial_out_file = $resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest.fa";
		my $cmd = "twoBitToFa -bed=$intialFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $initial_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive."BaseLine_one_closest_DHS/ClustersDHS/".$i."_DHS_closest_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $initial_out_file > $single";
 		system($cmd);
		
 	}	 	
} 

sub getUnlabeled_closestDHSperCluster
{
	$codeVersion = "Initialization_oneClosestGene/";
	my $mapfile = $dataDrive."Genome/FlyBase_Genes_map.tab";

	my %flymap;
	open IN, "<$mapfile" or die "Can not open file :$mapfile";

	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $annotationSymbo, $symbol) = split(/\s+/, $line);
 		push (@{$flymap{$ensemble}}, $flybaseID);
 	
 	}	
 	close IN;
	
	#get DHS used in initialization
	my %DHS_Initial;
	for (my $j=1; $j<= 39; $j++)
	{
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/".$j."_DHS_closest_in_uniq.bed ";
		open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
		
			chomp($line);
			my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
			my ($gene,$chr,$id,$stageW,$stage,$state ) = split(/\_/,$dhs_id);
			my $theID = $chr."_".$id."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
		close IN;	
	}
	my $k = keys %DHS_Initial;
	print ("number of DHS in initialization = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_10bp_uniq.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/DHS_unlabeled.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if( exists $DHS_Initial{$dhs_id})
		{
			print ("DHS $dhs_id exists\n");
		}
		else
		{
			print OUT  ("$line\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/DHS_unlabeled_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq_average";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=1;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$score,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			print OUTL ("$chr_dhs\t$s\t$e\t$flymap{$gene}[0]_".$dhs_id."_unlabeled\n");
			$index++;
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}

sub concatenate_ubiquitous_closestDHSperCluster
{
	$codeVersion = "BaseLine_one_closest_DHS/";
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_DHS_closest_in_uniq.bed ";	
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);

	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = "rm ". $resultDrive.$codeVersion."/ClustersDHS/".$j."_DHS_closest_single.fa ";
		system($cmd);
	
	}
	
}

#####Baseline_3 all DHS in +/- 50
sub all_DHS_50kb
{
	
	my $midDHSFile =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed";
	my $original =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";

	my %dhs;
	 
	open IN, "<$original" or die "Can't open file:$original";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	my $codeVersion= "BaseLine_all_DHS_50kb/";
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window_uniq";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		#print $cmd;
		system($cmd);
		
	}
	for (my $i=1;$i<=39; $i++)
	{
		
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window_non_uniq";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		#print $cmd;
		system($cmd);
	
		$cmd = "cat ".$resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window_uniq ".$resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window_non_uniq > ".$resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window";
		system($cmd);
		
	}

	for (my $i=1;$i<=39; $i++)
	{
		
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_50kb_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$genegc,$geneflybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			
			print OUTL ("$chr_dhs\t$s\t$e\t$geneflybase"."_".$dhs_id."_50kb\n");
			
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_50kb_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
	
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_50kb.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_50kb_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_50kb.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
		
	}	
	
	print ("\n");
 	
	
}

sub getUnlabeled_all_50kb
{
	$codeVersion = "BaseLine_all_DHS_50kb/";
	my $mapfile = $dataDrive."Genome/FlyBase_Genes_map.tab";

	my %flymap;
	open IN, "<$mapfile" or die "Can not open file :$mapfile";

	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $annotationSymbo, $symbol) = split(/\s+/, $line);
 		push (@{$flymap{$ensemble}}, $flybaseID);
 	
 	}	
 	close IN;
	
	#get DHS used in initialization
	my %DHS_Initial;
	for (my $j=1; $j<= 39; $j++)
	{
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/".$j."_DHS_closest_in_uniq.bed ";
		open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
		
			chomp($line);
			my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
			my ($gene,$chr,$id,$stageW,$stage,$state ) = split(/\_/,$dhs_id);
			my $theID = $chr."_".$id."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
		close IN;	
	}
	my $k = keys %DHS_Initial;
	print ("number of DHS in initialization = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_10bp_uniq.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/DHS_unlabeled.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if( exists $DHS_Initial{$dhs_id})
		{
			print ("DHS $dhs_id exists\n");
		}
		else
		{
			print OUT  ("$line\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/DHS_unlabeled_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq_average";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=1;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$score,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			print OUTL ("$chr_dhs\t$s\t$e\t$flymap{$gene}[0]_".$dhs_id."_unlabeled\n");
			$index++;
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}

sub concatenate_ubiquitous_all_50kb
{
	$codeVersion = "BaseLine_all_DHS_50kb/";
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_50kb_in_uniq.bed ";	
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);

	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = "rm ".$resultDrive.$codeVersion."/ClustersDHS/".$j."_50kb_single.fa ";
		system($cmd);	
	}
	
}
########################Second stage MC
sub getTSSperCluster_second
{
	for (my $i=1; $i<=39;$i++)
 	{
 		my $cluster_file = $dataDrive . "/geneClusters/$i.gene.cluster.non.uniq";		
 		open IN, "<$cluster_file" or die "Can not open file :$cluster_file";
 		my %genes;
 		my $line;
 		while ($line = <IN> )
	 	{
	 		chomp($line);
	 		my ($cg_gene, $flybase) = split(/\s+/, $line);
	 		$genes{$flybase} = 1;
	 	}
	 	close IN;
	 	my $noOfGenes = keys %genes;
 		
		my $tssFile = $dataDrive. "Genome/Gene_TSSs_2bp.bed";
		open IN, "<$tssFile" or die "Can not open file :$tssFile";
		my $tss_geneFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		open OUT, ">$tss_geneFile" or die "Can't open file :$tss_geneFile";
		my $line;
		my $count=0;
		my @remove;
	 	while ($line = <IN> )
	 	{
	 		chomp($line);
	 		
	 		my ($chr,$start,$end,$gene,$flybaseID,$strand) = split(/\s+/, $line);
	 		if (exists $genes{$flybaseID})
	 		{
	 			print OUT ("$line\n");
	 			$count++;
	 			push (@remove, $flybaseID);
	 			
	 		}
	 		
	 	}
	 	my %removed;
		@removed{ @remove } = delete @genes{ @remove };
	
		my $leftKeys = keys %genes;
		
	 	if($leftKeys ==0)
	 	{
	 		print ("cluster $i is done nicely\n");
	 	}
	 	else
	 	{
	 		print ("cluster $i genes = $noOfGenes, found only is $count, still remaining$leftKeys\n");
	 		print Dumper \%genes;
	 	}
	 	close OUT;
 	}
}

sub filterClusterTSS
{
	#for each gene, only write down one TSS which is the average of all the TSS's reported for this gene
	for (my $i=1; $i<=39;$i++)
 	{
		my $tss_geneFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		open IN, "<$tss_geneFile" or die "Can't open file :$tss_geneFile";
		
		my $line;
		my %gene_info;
		my $prev_gene="";
		while ($line = <IN> )
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$v,$strand) = split(/\s+/, $line);
			$gene_info{$gene}{'chr'} = $chr;
			$gene_info{$gene}{'strand'}=$strand;
			$gene_info{$gene}{'start'} +=$start;
			$gene_info{$gene}{'end'} +=$end;
			$gene_info{$gene}{'count'}++;
		}
		close IN;
		
		my $tss_median_geneFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq_average";
		open OUT, ">$tss_median_geneFile" or die "Can't open file :$tss_median_geneFile";
		
		foreach my $gene (keys %gene_info)
		{
			if($gene_info{$gene}{'strand'} eq "+")
			{
				my $start = int ($gene_info{$gene}{'start'} / $gene_info{$gene}{'count'});
				my $end = $start +1;
				print OUT ("$gene_info{$gene}{'chr'}\t$start\t$end\t$gene\t100\t$gene_info{$gene}{'strand'}\n");
			}
			else
			{
				my $end = int ($gene_info{$gene}{'end'} / $gene_info{$gene}{'count'});
				my $start = $end -1;
				print OUT ("$gene_info{$gene}{'chr'}\t$start\t$end\t$gene\t100\t$gene_info{$gene}{'strand'}\n");
			}
		}
		close OUT;
		
 	}
}

sub splitDHSIntoClusters_specifiedGenes_2nd_stage
{
	mapGeneIDs();
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id,@rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_specified_genes.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		 my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$gene) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	#my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_specified_genes_used_in_initialization.bed";
	#open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=10; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			
 			my $gene = $flymap{$line};
 		
 			if(exists($gene_DHS{$gene}))
	 		{
	 				
	 			for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 			{
	 				my $s = $gene_start{$gene}[$j];
	 				my $e = $gene_end{$gene}[$j] ;
	 				print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 			#	print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 				$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 			}
	 			$geneCount{$i}{$gene}++;
	 		}
	 			
 		}
 	
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $clusterFile_fa > $uniq_file";
 		#system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		#system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		#system($cmd);
 	
 	}
 	close ALL;
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
		
 
}

sub getClosestPerClusterForUnspecifiedGenes_2nd_stage
{
	#intersectBed -a DHS_no_overlap_TSS_10bp.bed -b ../../CRM/redfly_all_unspecified.bed -u  > DHS_overlap_redfly_unspecified_genes.bed
	#closestBed -a DHS_overlap_redfly_unspecified_genes_midpoints.bed -b ../../Genome/fly_TSS_10bp_uniq.bed -d -t all > DHS_overlap_redfly_unspecified_genes_closest.bed
	#read mapping data
	#flybase ID could be linked to more than one gene
	my $mapfile = $dataDrive."Genome/FlyBase_Genes_map.tab";

	my %flymap;
	open IN, "<$mapfile" or die "Can not open file :$mapfile";

	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $annotationSymbo, $symbol) = split(/\s+/, $line);
 	
 
 		push (@{$flymap{$ensemble}}, $flybaseID);
 	
 	}	
 	close IN;
	
	
	
	my %geneClusters;
	my $line;
	for (my $i=1;$i<=10; $i++)
	{
		my $clusterFile = $dataDrive."geneClusters/".$i.".gene.cluster";
		open IN, "<$clusterFile" or die "Can't open file:$clusterFile";
		while ($line = <IN>)
		{
			chomp($line);
			#my @flybase = $flymap{$line};
			push (@{$geneClusters{$i}}, $line);
		}
	}
	close IN;
	
	#to get the actual start and end
	my %dhs;
	my $dhsFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes.bed";
	 
	open IN, "<$dhsFile" or die "Can't open file:$dhsFile";
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	my $closestFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes_closest.bed";
	open IN, "<$closestFile" or die "Can't open file:$closestFile";
	#my %geneClusters_clone = dclone(%geneClusters);
	my %closestGenesDHS;
	my %GeneCount;
	my %DHSCount;
	my %dhs_gene;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$dhs_mid_start,$dhs_mid_end,$dhs_id,$chr_g,$gene_tss_s,$gene_tss_e,$gene,$score, $strand,$distance) = split(/\s+/, $line);
		
		for (my $i=1;$i<=10; $i++)
		{	
			#my $flybase = $flymap{$gene}[0];		
			my @output = grep {$gene eq $_} @{$geneClusters{$i}};
			if(@output >0)
			{
				$DHSCount{$i}{$dhs_id}=1;
				$GeneCount{$i}{$gene}++;
				push(@{$closestGenesDHS{$i}}, $dhs_id);
				$dhs_gene{$dhs_id}=$gene;
				#print ("$gene \t cluster $i\n");
			}
		}
		
	}
	
	for (my $i=1;$i<=10; $i++)
	{
		my $c = keys %{$GeneCount{$i}};
		my $s = keys %{$DHSCount{$i}};
		print ("cluster $i\t$c\t$s\n");
	}
	

	my $allUNSPEC = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/redfly_dhs_unspecified_genes_used_in_initialization.bed";
	open UNSPEC, ">$allUNSPEC" or die "Can't open file:$allUNSPEC";
	
	for (my $i=1;$i<=10; $i++)
	{
		my $index=1;
		my $outputFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_closest_dhs_unspecified_genes.bed";
		open OUT, ">$outputFile" or die "Can't open file:$outputFile";
		if(exists $closestGenesDHS{$i})
		{
			for(my $l=0; $l <@{$closestGenesDHS{$i}}; $l++)
			{
				my $id = $closestGenesDHS{$i}[$l];
				my $s = $dhs{$id}{'start'};
				my $e = $dhs{$id}{'end'};
				print OUT ("$dhs{$id}{'chr'}\t$s\t$e\t$flymap{$dhs_gene{$id}}[0]_".$id."_initial\n");
				print UNSPEC ("$dhs{$id}{'chr'}\t$dhs{$id}{'start'}\t$dhs{$id}{'end'}\t$id\tcluster_$i\t$flymap{$dhs_gene{$id}}[0]\n");
				$index++;
			}
		}
		close OUT;
	
		my $redFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly.bed";
		my $allFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_initialize_redfly.bed";
		my $cmd = "cat $redFile $outputFile > $allFile";
		system($cmd);
		
		my $mergedFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_initialize_redfly_uniq.bed";
		$cmd= "sort -k1,1 -k2,2n -k3,3n --uniq $allFile > $mergedFile";	
		system($cmd);
		

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_initialize_redfly_uniq.fa";
		my $cmd = "twoBitToFa -bed=$mergedFile -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_initialize_redfly_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
	}
	
		close UNSPEC;
		
		my $cat = "cat ".$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/redfly_dhs_specified_genes_used_in_initialization.bed ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_unspecified_genes_used_in_initialization.bed > ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_all_used_in_initialization.bed";
		system($cat);
		
		
}

sub concatenate_ubiquitous_2nd_stage
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$j."_specified_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

sub getUnlabeled_2nd_stage
{
	#get DHS used for uniq genes
	my %DHS_Initial;
	for (my $j=11; $j<= 39; $j++)
	{
		
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$j."_pos.fa";
		open IN, "<$initializeFile" or print "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
			chomp($line);
			if($line =~ "^\>")
			{
			
				my ($gene,$chr,$start,$stageW, $stage,$status) = split(/\_/, $line);
				my $theID = $chr."_".$start."_".$stageW."_".$stage;
				$DHS_Initial{$theID}++;
			}
		}
		close IN;	
	}
	
	#read ubiquitous
	my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all_single.fa";
	open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		if($line =~ "^\>")
		{
			my ($gene,$chr,$start,$stageW, $stage,$status) = split(/\_/, $line);
			my $theID = $chr."_".$start."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
	}
	close IN;	
	my $k = keys %DHS_Initial;
	print ("number of DHS after MC first stage = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/DHS_unlabeled_second_stage.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if(exists $DHS_Initial{$dhs_id})
		{
			#print ("DHS $dhs_id does  exist\n");
		}
		else{
			print OUT  ("$line\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/DHS_unlabeled_second_stage_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=11;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=11;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$flybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			print OUTL ("$chr_dhs\t$s\t$e\t$flybase"."_".$dhs_id."_unlabeled2\n");
			$index++;
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}

sub splitDHSIntoClusters_TSS_all_genes
{
	mapGeneIDs();
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$cg_gene, $gene, $strand) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_TSS.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			
 			my $flybase = $flymap{$line};
 				if(exists($gene_DHS{$flybase}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
	 				{
	 					my $s = $gene_start{$flybase}[$j];
	 					my $e = $gene_end{$flybase}[$j] ;
	 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$flybase"."_".$gene_DHS{$flybase}[$j]."_TSS\n");
#	 					print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$flybase}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$flybase}++;
	 			
 			}
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_TSS_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_TSS_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_TSS_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
}

sub concatenate_ubiquitous_TSS_all_genes
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_TSS.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$j."_TSS_uniq.bed ";  # _initialize_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_TSS_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_TSS.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_TSS_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}
