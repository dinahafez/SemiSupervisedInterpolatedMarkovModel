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
my %mapping;
mapGeneIDs();
writeOneFlyBaseID();
writeRedFly();
my %genes;
my %genes_flybase;

sub mapGeneIDs
{
	my $mapfile = $dataDrive."Genome/Gene_map_new.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr, $start, $end, $strand) = split(/\s+/, $line);
 		if(! grep $_ eq $flybaseID, @{$flymap{$ensemble}})
 		{
 			push (@{$flymap{$ensemble}}, $flybaseID);
 		}
 	}	
 	close IN;
	foreach my $gene (keys %flymap)
	{
		if ($flymap{$gene} > 1)
		{
			for(my $k=1; $k< @{$flymap{$gene}}; $k++)
			{
				$mapping{$flymap{$gene}[$k]} =$flymap{$gene}[0]; 
			}
		}
	}
}

sub writeOneFlyBaseID
{
	my $mapfile = $dataDrive."Genome/Gene_map_new.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $outFile = $dataDrive."Genome/Gene_map_new_filtered.tab";
	open OUT, ">$outFile" or die "Can not open file :$outFile";
	
	my $line = <IN> ;
	my $count=0;
	my %bad;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr, $start, $end, $strand) = split(/\s+/, $line);
 	
 		if ($ensemble ne "CG11861")
 		{
	 		if (@{$flymap{$ensemble}} > 1)
	 		{
	 			print OUT ("$ensemble\t${$flymap{$ensemble}}[0]\t$chr\t$start\t$end\t$strand\n");
	 			if(exists $genes_flybase{$flybaseID} && ($genes_flybase{$flybaseID} ne $flymap{$ensemble}[0]))
	 			{
	 				print ("$flybaseID\t$genes_flybase{$flybaseID}\t$flymap{$ensemble}[0]\n");
	 			}
	 			else
	 			{
	 				
	 				$genes_flybase{$flybaseID} = $flymap{$ensemble}[0];
	 			}
	 			
	 			$bad{$ensemble} =1;			
	 		}
	 		else
	 		{
	 			print OUT ("$line\n");
	 		}
 		}
 	}
 		
 	close IN;
	
	my $k = keys %bad;
 	print ("no of genes that have multiple flybase (only 2) = $k\n");
 
 my $cmd = "sort -k1,1  -k2,2 -k3,3 -k4,4n -k5,5 -k6,6 Gene_map_new_filtered.tab --uniq > Gene_map_new_filtered_uniq.tab";
 #system ($cmd);

}

sub writeRedFly
{
	#write redfly using the same FBgn ID
	my $inFile = $dataDrive."CRM/redfly_all_specified.bed";
	open IN, "<$inFile" or die "Can not open file :$inFile";
	
 	my $outFile = $dataDrive."CRM/redfly_all_specified_modified.bed";
	open OUT, ">$outFile" or die "Can not open file :$outFile";
	
 	while (my $line = <IN> )
 	{
 		chomp($line);
 		my ($chr, $start, $end, $gene) = split(/\s+/, $line);
 		if (exists $genes_flybase{$gene})
 		{
 			print OUT ("$chr\t$start\t$end\t$genes_flybase{$gene}\n");
 		}
 		else
 		{
 			print OUT ("$line\n");
 		}
 	
 	}
 	close IN;
 	close OUT;
}

sub mapGeneIDs2
{
	my $mapfile = $dataDrive."Genome/Gene_map_new.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $line = <IN> ;
	my %bad;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr,$start, $end, $strand) = split(/\s+/, $line);
 		if(exists $flymap{$ensemble} && ($flymap{$ensemble} ne $flybaseID))
 		{
 			if(exists $bad{$ensemble} && ($bad{$ensemble} ne $flybaseID))
 			{
 				print ("$ensemble\t$flymap{$ensemble}\t$bad{$ensemble}\t$flybaseID\n");  #does not happen
 			}
 			else
 			{
 				$bad{$ensemble} = $flybaseID;
 				#print ("$ensemble\t$flymap{$ensemble}\t$bad{$ensemble}\t$flybaseID\n");
 			}
  		}
 		else
 		{			
 				$flymap{$ensemble} = $flybaseID;
 		}
 	}	
 	
 	my $k = keys %bad;
 	print ("no of genes that have multiple flybase (only 2) = $k\n");
 	close IN; 
 	
 	my %chr;
 	my %start;
 	my %end;
 	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $line = <IN> ;
	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr, $start, $end, $strand) = split(/\s+/, $line);
 		if(exists $bad{$ensemble})
 		{
 			if(exists $chr{$ensemble} && ($chr{$ensemble} ne $chr))
 			{
 				print ("$ensemble\t$chr{$ensemble}\t$chr\n");
 			}
 			elsif(exists $start{$ensemble} && (abs($start{$ensemble} - $start)> 1000))
 			{
 				print ("$ensemble\t$start{$ensemble}\t$start\n");
 			}
 			
 			else
 			{
 				$chr{$ensemble} = $chr;
 				$start{$ensemble} = $start;
 				$end{$ensemble} = $end;
 			}
 		}
 	}
 		
}

