#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($file_in, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses S2N's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to S2N results file *
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


my ($counter, @temp, $gene, $pval, $score);
my $pval_threshold = 0.05;
$counter = 1;
open(IN, "tail -n+2 $file_in | cut -f 1,4,5 | grep -vw \"NA\" | sort -k3,3g -k2,2gr |");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tpValue\tScore\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[1];
	$pval = $temp[2];
	last if($pval > $pval_threshold);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $pval . "\t" . $score . "\n";
	$counter++;
}
close(OUT);
close(IN);
