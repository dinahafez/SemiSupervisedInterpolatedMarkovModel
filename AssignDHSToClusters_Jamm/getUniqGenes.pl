use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/";

#my $dataDrive = "../../data/";
#my $resultDrive = "../../results/";

my %flymap;
#mapGeneIDs();

my %genes;
my %genes_flybase;

for (my $i=1;$i<=39; $i++)
{
	my $geneFile = $dataDrive."geneClusters/".$i.".gene.cluster";
	open IN, "<$geneFile" or die "Can't open file:$geneFile";
 	my $line;
 
 	while ($line = <IN>)
	{
		chomp($line);
		#$genes{$line}++;
		my $flybaseID = $flymap{$line};
		$genes_flybase{$flybaseID}{$i}++;
		
		if($flybaseID eq "FBgn0041094")
		{
			print ("cluster is $i\n");
		}
	}
	close IN;
}
my $uniq_broad=0;
my $uniq_rest=0;
my $uniq=0;

for (my $i=1;$i<=0; $i++)
{
	my $geneFile = $dataDrive."geneClusters/".$i.".gene.cluster";
	my $geneFileUniq = $dataDrive."geneClusters/".$i.".gene.cluster.uniq";
	my $geneFileNonUniq = $dataDrive."geneClusters/".$i.".gene.cluster.non.uniq";
#	open IN, "<$geneFile" or die "Can't open file:$geneFile";
#	open OUT, ">$geneFileUniq" or die "Can't open file:$geneFileUniq";
#	open OUT2, ">$geneFileNonUniq" or die "Can't open file:$geneFileNonUniq";
 	my $line;
 
 	while ($line = <IN>)
	{
		chomp($line);
		my $flybaseID = $flymap{$line};
		my $noOfClusters = keys %{$genes_flybase{$flybaseID}};
		
		if ($line ne "CG11861")
		{
			if ($noOfClusters == 1)
			{
				print OUT ("$line\t$flybaseID\n");
				
			}
			else
			{
				print OUT2 ("$line\t$flybaseID\n");
			}
		}
		
	}
	close IN;
	close OUT;
	close OUT2;
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