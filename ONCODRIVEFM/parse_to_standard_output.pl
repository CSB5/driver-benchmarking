#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($file_in, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses OncodriveFM's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to OncodriveFM results file *
	--out = path to output file *
	--debug: prints trace to STDERR
	--help : prints this message 
	
* indicates required parameters	


Version:
	1.0

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}

GetOptions(
	"in=s"      	=> \$file_in,
	"out=s"         => \$file_out,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "INPUT: $file_in\n";
	print STDERR "OUTPUT: $file_out\n";	
}


my ($counter, @temp, $gene, $qval);
my $qval_threshold = 0.05;
$counter = 1;
open(IN, "sort -k2,2g $file_in |");
open(OUT, "> $file_out");
<IN>;	# skip header
<IN>;	# skip header
<IN>;	# skip header
<IN>;	# skip header
<IN>;	# skip header
print OUT "Gene_name\tSample\tRank\tqValue\tInfo\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$qval = $temp[2];
	last if($qval > $qval_threshold);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $qval . "\t" . "-" . "\n";
	$counter++;
}
close(OUT);
close(IN);
