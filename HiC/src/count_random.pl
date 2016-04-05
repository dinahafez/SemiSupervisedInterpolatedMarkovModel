#read all DHSs
#read genes TSS 

use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/";
my $codeVersion = "AssignDHSToRandomGene/permute/";

for (my $i = 1; $i<2; $i++)
{
	RandomSample($i);
	runCommands($i);
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

	my $windowFile = $resultDrive. $codeVersion."DHStoAll50kb_uniq_gene_2.bed";
	open IN, "<$windowFile" or die "Can't open file:$windowFile";

	my $outFile = $resultDrive. $codeVersion."DHStoRandomGene_$i.bed";
	open OUT, ">$outFile" or die "Can't open file:$outFile";



	my $prevId = "";
	my @genes=();
	my $sum =0 ;
	my $total =0;

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
				my $noOfGenes = @genes;
				$sum = $sum + $noOfGenes;
				$total = $total +1;
				#my $geneNo = int(rand (@genes));
				#print OUT ("$dhs{$dhs_id}{'chr'}\t$dhs{$dhs_id}{'start'}\t$dhs{$dhs_id}{'end'}\t$genes[$geneNo]\n");
			}
			#do cleanups
			@genes = ();
			push (@genes, $geneFB);
			$prevId = $dhs_id;
		}

	}
	print ("Total number of DHSs = $total\n);
	
close IN;
close OUT;
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


