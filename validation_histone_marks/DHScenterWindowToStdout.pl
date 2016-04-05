#!/usr/bin/perl

use strict;
use warnings;


use feature ':5.10';
use Data::Dumper;
use List::MoreUtils qw/ uniq /;



my @commands = @ARGV;
use Cwd;
my $dir = getcwd;

open INPUT, "<$commands[0]" or die "Cannot open input file.\n";
#open OUT, ">$dir/$commands[0].strandedDHSwindow" or die "Cannot open output file.\n";
my %sites;


while (<INPUT>) 
{
	my @d = split;
	my $center = sprintf("%.0f", ($d[2] - (($d[2] - $d[1])/2))); 
	my $start = $center - $commands[1];
	next if $start < 0;
	my $end = $center + ($commands[2] + 1);
	my $size = $d[2] - $d[1];
	print "$d[0]\t$start\t$end\t$size\n"; #\t$d[3]\n";
}
