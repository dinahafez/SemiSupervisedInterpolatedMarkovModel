 #!/usr/bin/perl

use strict;
use warnings;


use feature ':5.10';
use Data::Dumper;
use List::MoreUtils qw/ uniq /;


use Cwd;
my $dir = getcwd;
my @commands = @ARGV;


open IP, "<$commands[0]" or die "Cannot open input file.\n";
open WINDOW, "<$commands[1]" or die "Cannot open input file 2.\n";
open OUT, ">$commands[2].peakCoverageMeta" or die "Cannot open output file.\n";

my %tss;
my %input;
my %ip;
my @tss;

my $binsize = $commands[3];

while (<IP>)
{
	chomp;
	my @c = split;
	for (my $i = $c[1]; $i < $c[2]; $i++)
	{
		$ip{"$c[0]\t$i"} = 1;
	}		
}


my $id_count=0;

while (<WINDOW>) #for every bp position in window create strand specific hash of arrays where keys are ClusterID_Clustershape_ClusterRelativeSense and array is length of window with score, if exists, or "0" at each position  
{
	chomp;
	my @d = split;
		my $bincount = 0;
		my $ipbinscore = 0;
$id_count++;
		#push @{$tss{"$d[6]_$d[3]"}}, $d[7];
		for (my$i=$d[1]; $i<$d[2]; $i++)
		{
			$bincount++;
			if (exists $ip{"$d[0]\t$i"}) 
			{
				$ipbinscore = $ip{"$d[0]\t$i"};
			}		
			if ($bincount == $binsize)
			{
				push @{$tss{"dhs_".$id_count}}, $ipbinscore;
				$bincount = 0;
				$ipbinscore = 0;
				
			} 
}
}


foreach my $keys (sort keys %tss)
{
	print OUT $keys, "\t";
	for (my$i=0; $i < @{$tss{$keys}}; $i++)
	{
		print OUT $tss{$keys}[$i], "\t";
	}
	print OUT "\n";
}


