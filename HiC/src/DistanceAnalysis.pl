use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);
use List::MoreUtils qw/pairwise/;

my $i = $ARGV[0];
my $j=$ARGV[1];

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/"; 

#read all TSS
my %Tss;
readTSS ();

my $originalFile = $resultDrive."model_intersect_partnerTogene.bed.filtered.uniq.bed";
my $overlapFile=$resultDrive."model_overlap_HiC_uniq.bed";
my $overlapDist=$resultDrive."model_overlap.dist";
my $nonOverlapDist=$resultDrive."model_no_overlap.dist";

my %overlap;	
#read overlap
open IN, "<$overlapFile" or die "Can not open file :$overlapFile";
 
while (my $line = <IN> )
{
	
	chomp($line);
	my @fields = split(/\s+/,$line);
	push(@{$overlap{$fields[3]}}, $fields[1]);
}
close IN;


open IN, "<$originalFile" or die "Can not open file :$originalFile";
open OUT, ">$overlapDist" or die "Can not open file :$overlapDist";
 open OUT2, ">$nonOverlapDist" or die "Can not open file :$nonOverlapDist";
while (my $line = <IN> )
{
		
	chomp($line);
	
	my @fields = split(/\s+/,$line);
	
	#get minumum distance between DHS and TSS
	my $mid = int((($fields[2] - $fields[1]) / 2)	+ $fields[1]);
	my @mids  = ($mid) x @{$Tss{$fields[3]}};
	my @distances =abs( pairwise { $a -  $b } @{$Tss{$fields[3]}}, @mids);
	#my @distances =abs(@{$Tss{$fields[3]}}- @mids);
	my @sorted = sort {$a <=> $b} @distances; 
	my $dist = int($sorted[0]);

	#my $dist=10000000;
	#for (my $k=0; $k< $Tss{$fields[3]}; $k++)
	#{
	#	my $d = int(abs ($Tss{$fields[3]}[$k ] - $mid));
	#	if ($d < $dist )
	#	{
	#		$dist = $d;
	#	}
	#}	
	
	#print ("$dist\n");
	if (grep(/^$fields[1]/i, @{$overlap{$fields[3]}}))
	{#overlaps
		#print ("it overlaps \n");
		print OUT ("$line\t$dist\n");
	}
	else
	{#print ("it does not overlaps \n");
		print OUT2 ("$line\t$dist\n");
	}

	if ($fields[3] eq "FBgn0031174")
	{print ("$line \t Distance $dist\n");

	print ("mids\n");
	print Dumper @mids;
	print ("distances\n");	
	print Dumper @distances;	
		}
}
close IN;
close OUT;
close OUT2;


sub readTSS
{
	print ("Reading TSSs\n");
	for (my $i=11;$i<=39; $i++)
	{
		my $bed = $dataDrive ."geneClusters/".$i."_cluster_TSS_uniq";
		open IN, "<$bed" or print "Can't open file:$bed\n";
	 	my $line;
	 
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene_name,$gene)= split(/\s+/, $line);
			#$Tss{$gene} = $start + 1;		
			push (@{$Tss{$gene}}, $start+1);
		}
		close IN;

		my $bed = $dataDrive ."geneClusters/".$i."_cluster_TSS_non_uniq";
		open IN, "<$bed" or print "Can't open file:$bed\n";	 	
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene_name,$gene,$strand)= split(/\s+/, $line);
			#$Tss{$gene} = $start + 1;
			
			push (@{$Tss{$gene}}, $start+1);
			if ($gene eq "FBgn0031174")
			{print ("$gene\t$start\n");}
		}
		
		close IN;
	}


}
	
 		
