#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($file_in, $file_out, $flag_help, $flag_debug) = @ARGV;

my $help_message = "
This script parses oncoIMPACT's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to DriverNet results file *
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


my ($counter, @temp, $gene, $score);
my $score_threshold = 0;
$counter = 1;
open(IN, "sort -k15,15nr $file_in |");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[14];
	last if($score <= $score_threshold);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
	$counter++;
}
close(OUT);
close(IN);
