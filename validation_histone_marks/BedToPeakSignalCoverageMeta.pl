#!/usr/bin/perl

use strict;
use warnings;


use feature ':5.10';
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

#This script takes as input one tab delimited "score" file as created by 5end_read_bed.pl consisting of single nucleotide resolution read counts as
#chr\start\end\strand\score and a second tab delimited "window" file as created by ClustsToWindow.pl consisting of windows as 
#chr\start\end\name\strand, assigns the score at every position from the "score" file to the appropriate window and position, 
#outputing sense and antisense nucleotide-resolution window scores for every name as tab delimited rows 
use Cwd;
my $dir = getcwd;
my @commands = @ARGV;


open IP, "<$commands[0]" or die "Cannot open input file.\n"; #reads 
open SUMMIT, "<$commands[1]" or die "Cannot open input file.\n"; #peaks
open WINDOW, "<$commands[2]" or die "Cannot open input file.\n"; #dhs_window file
print ("$commands[3].$commands[5]bin.peakSignalCoverageMeta\n");
open OUT, ">$commands[3].$commands[5]bin.peakSignalCoverageMeta" or die "Cannot open output file.\n";

my %tss;
my %input;
my %summit;
my %ip;
my @tss;

my $ipfragment = $commands[4];
my $binsize = $commands[5];


while (<IP>)
{
        chomp;
        my @c = split;
        if ($c[5] eq "+")
        {
                for (my $i = $c[1]; $i < $c[1] + $ipfragment; $i++)
                {
                        if (exists $ip{"$c[0]\t$i"})
                        {
                                $ip{"$c[0]\t$i"}++;
                        }
                        else
                        {
                                $ip{"$c[0]\t$i"} = 1;
                        }
                }
        }
        elsif ($c[5] eq "-")
        {
                for (my $i = $c[2] -1; $i >= $c[2] - 1 - $ipfragment; $i--)
                {
                        if (exists $ip{"$c[0]\t$i"})
                        {
                                $ip{"$c[0]\t$i"}++;
                        }
                        else
                        {
                                $ip{"$c[0]\t$i"} = 1;
                        }
                }
        }
}



while (<SUMMIT>)
{
	chomp;
	my @c = split;
	for (my $i = $c[1]; $i < $c[2]; $i++)
	{
		$summit{"$c[0]\t$i"} = 1;
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
			if ((exists $summit{"$d[0]\t$i"}) && (exists $ip{"$d[0]\t$i"}))
			{
				$ipbinscore += $ip{"$d[0]\t$i"};
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


