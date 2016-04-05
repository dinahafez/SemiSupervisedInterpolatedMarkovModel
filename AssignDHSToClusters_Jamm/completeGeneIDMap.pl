use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

#This script does the mapping for the gene names from CGI index to FBGN name 
my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";



my %flymap;
readMappedGenes();
getUnMappedGenes();
#doMap();

sub readMappedGenes
{
	#read the  map file
	my $file = $dataDrive."Genome/Gene_map_new.tab";
	
	open IN, "<$file" or die "Can not open file :$file";
	my $line = <IN> ;
	while ($line = <IN> )
	{
		chomp($line);
	 	my ($ensemble, $flybaseID, $chr,$start, $end, $strand) = split(/\s+/, $line);
	 	$flymap{$ensemble} = $flybaseID;
	}
	close IN;
}

sub getUnMappedGenes
{
	for (my $i=1; $i<=39;$i++)
 	{
 		my $cluster_file = $dataDrive . "/geneClusters/$i.gene.cluster";		
 		open IN, "<$cluster_file" or die "Can not open file :$cluster_file";
 		my %genes;
 		my $line;
 		while ($line = <IN> )
	 	{
	 		chomp($line);
	 		if(!exists $flymap{$line})
	 		{
	 			print ("$line\n");
	 		}
	 		
	 	}
	}
}

sub doMap
{
	
	#read the  map file
	my $file = $dataDrive."Genome/Gene_map_new.tab";
	
	open IN, "<$file" or die "Can not open file :$file";
	my $line = <IN> ;
	while ($line = <IN> )
	{
		chomp($line);
	 	my ($ensemble, $flybaseID, $chr,$start, $end, $strand) = split(/\s+/, $line);
	 	$flymap{$flybaseID}{"cg"} = $ensemble;
	 	push(@{$flymap{$flybaseID}{"chr"}}, $chr);
	 	push(@{$flymap{$flybaseID}{"start"}},  $start);
	 	push(@{$flymap{$flybaseID}{"end"}}, $end);
	 	push(@{$flymap{$flybaseID}{"strand"}}, $strand);
	}
	close IN;
	
	my $missing = $dataDrive."Genome/map";
	my $rest = $dataDrive."Genome/missing_mapped_genes.txt";
	open IN, "<$missing" or die "Can not open file :$missing";
	open OUT, ">$rest" or die "Can not open file :$rest";
	my %wrote;
	my %miss;
	while ($line = <IN> )
	{
		chomp($line);
	 	my ($cgname, $flybase, @rest) = split(/\s+/, $line);
		if (exists $flymap{$flybase})
		{
			for(my $i=0; $i<@{$flymap{$flybase}{"chr"}}; $i++)
			{
				print OUT ("$cgname\t$flybase\t$flymap{$flybase}{\"chr\"}[$i]\t$flymap{$flybase}{\"start\"}[$i]\t$flymap{$flybase}{\"end\"}[$i]\t$flymap{$flybase}{\"strand\"}[$i]\n");
			}
			$wrote{"$cgname"}=1;
		} 
		else
		{
			#print ("$cgname\t$flybase\n");
			$miss{"$cgname"} = $flybase;
		}
	}
	
	for my $key (keys %miss)
	{
		if (!exists $wrote{$key})
		{
			print ("$key\t$miss{$key}\n");
		}
	}
	close IN; 
	close OUT;
}
