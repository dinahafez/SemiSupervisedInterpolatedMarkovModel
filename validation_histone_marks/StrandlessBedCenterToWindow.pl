#!/usr/bin/perl

use strict;
use warnings;


use feature ':5.10';
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

#This script takes three input commands. The first is a bed format file consisting of a position in the score field from which a window is to be made
#from that position minus command two to that position plus command three. the output is a bed format file of the produced window maintaining the name field
#from input file and replacing the score field with a counter name.


my @commands = @ARGV;
use Cwd;
my $dir = getcwd;

open INPUT, "<$commands[0]" or die "Cannot open input file.\n";
open OUT, ">$dir/$commands[0]_$commands[1]up_$commands[2]down.bedCenterWindow" or die "Cannot open output file.\n";

my $count = 0;
while (<INPUT>) 
{
	$count++;
	chomp;
	my @d = split;
	my $length = $d[2] - $d[1];
	my $center = sprintf("%.0f", ($d[1] + (($d[2] - $d[1])/2))); 
	my $start = $center - $commands[1];
	next if $start < 0;
	my $end = $center + ($commands[2] + 1);
	print OUT "$d[0]\t$start\t$end\t$count\t\.\t\+\t$count\n";
}
