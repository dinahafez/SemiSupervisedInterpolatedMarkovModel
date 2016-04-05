#read all DHSs
#read genes TSS 

use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/";
my $codeVersion = "permute/";

for (my $i = 200; $i<=300; $i++)
{
	RandomSample_FromGene($i);
	#runCommands($i);
	
}

	
sub RandomSample
{
	my ($i) = @_;

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
	close IN;

	my $windowFile = $resultDrive. $codeVersion. "DHStoAll50kb_uniq_gene_2_sorted.bed";  #$codeVersion."DHStoAll50kb_uniq_gene_2.bed";
	open IN, "<$windowFile" or die "Can't open file:$windowFile";

	my $outFile = $resultDrive. $codeVersion."DHStoRandomGene_$i.bed";
	open OUT, ">$outFile" or die "Can't open file:$outFile";



	my $prevId = "";
	my @genes=();
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id,$chrg,$startg,$endg, $geneSymbol, $geneFB, $strand) = split(/\s+/, $line);
		if ($dhs_id eq $prevId)
		{
			push (@genes, $geneFB);
		}
		else
		{
			#randomly select a gene
			if ($prevId ne "")
			{
				my $geneNo = int(rand (@genes));
				print OUT ("$dhs{$dhs_id}{'chr'}\t$dhs{$dhs_id}{'start'}\t$dhs{$dhs_id}{'end'}\t$genes[$geneNo]\n");
			}
			#do cleanups
			@genes = ();
			push (@genes, $geneFB);
			$prevId = $dhs_id;
		}

	}
close IN;
close OUT;
	my $cmd = "./getOverlapHiC_FromBed.sh ".$codeVersion."DHStoRandomGene_".$i.".bed 2";
	print ("$cmd\n");
	system($cmd);
}

# windowBed -a genes_in_analysis.uniq.non.uniq -b /data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed -w 50000 > geneTo50kDHS2.bed	
sub RandomSample_FromGene
{
	my ($i) = @_;

	my $midDHSFile =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed";
	my $original =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";

	my %DHS;
	 
	open IN, "<$original" or die "Can't open file:$original";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$DHS{$dhs_id}{'start'} = $start;
		$DHS{$dhs_id}{'end'} = $end;
		$DHS{$dhs_id}{'chr'} = $chr;
		
	}
	close IN;

	my $windowFile = $resultDrive. "geneTo50kDHS_sorted2.bed ";  
	open IN, "<$windowFile" or die "Can't open file:$windowFile";

	my $outFile = $resultDrive. $codeVersion."DHStoRandomGene_$i.bed";
	open OUT, ">$outFile" or die "Can't open file:$outFile";



	my $prevId = "";
	my @dhs=();
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$geneSymbol, $geneFB, $strand,$chrg,$startg,$endg, $dhs_id) = split(/\s+/, $line);
		if ($geneFB eq $prevId)
		{
			push (@dhs, $dhs_id);
			
		}
		else
		{
			#randomly select a gene
			if ($prevId ne "")
			{
				for (my $j=0; $j<6; $j++)
				{
					my $no = int(rand (@dhs));
					my $dhs_id = $dhs[$no];
					print OUT ("$DHS{$dhs_id}{'chr'}\t$DHS{$dhs_id}{'start'}\t$DHS{$dhs_id}{'end'}\t$geneFB\n");
				}
			}
			#do cleanups
			@dhs = ();
			push (@dhs, $dhs_id);
			$prevId = $geneFB;
		}

	}
close IN;
close OUT;
	my $cmd = "./getOverlapHiC_FromBed.sh ".$codeVersion."DHStoRandomGene_".$i.".bed 2";
	print ("$cmd\n");
	system($cmd);
}

sub RandomSample_AnyGene
{
	my ($i) = @_;

	my $midDHSFile =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed";
	my $original =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";

	my %dhs;
	my %genes;
	my @gene_array;

	my $geneFile = $resultDrive."all_clusters_TSS_uniq_and_non_uniq";
	open IN, "<$geneFile" or die "Can't open file:$geneFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$symbol, $gene_id, $strand) = split(/\s+/, $line);
		$genes{$gene_id}{'start'} = $start;
		$genes{$gene_id}{'end'} = $end;
		$genes{$gene_id}{'chr'} = $chr;
		push (@gene_array,$gene_id);
		
	}
	close IN;
	
	my $outFile = $resultDrive. $codeVersion."DHStoRandomGene_$i.bed";
	open OUT, ">$outFile" or die "Can't open file:$outFile";
	open IN, "<$original" or die "Can't open file:$original";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		my $geneNo = int(rand (@gene_array));
		my $gene = $gene_array[$geneNo];
		print OUT ("$chr\t$start\t$end\t$gene\n");
		
		
	}
	close IN;


close OUT;
	my $cmd = "./getOverlapHiC_FromBed.sh ".$codeVersion."DHStoRandomGene_".$i.".bed 2";
	print ("$cmd\n");
	system($cmd);
}

sub runCommands
{
	my ($i) = @_;
	my $outFile = $resultDrive. $codeVersion."DHStoRandomGene_$i.bed";

	my $cmd = "intersectBed -a $outFile -b ".$dataDrive."HiC/sextonCell_95_partnerTogene.bed -u  > ".$outFile."_intersect_sexton_partnerTogene.bed";
	system ($cmd);

	$cmd = "perl fiterInteractions.pl ".$outFile."_intersect_sexton_partnerTogene.bed";
	system ($cmd);
	$cmd = "rm $outFile";
	system ($cmd);
	$cmd = "rm ".$outFile."_intersect_sexton_partnerTogene.bed";
	system ($cmd);


}


