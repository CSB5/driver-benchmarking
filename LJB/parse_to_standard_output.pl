#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($file_in, $outDir, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses ANNOVAR's annotation output to a standard output for LJB methods.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to ANNOVAR's annotation results file *
	--outDir = path to output folder *
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
	"outDir=s"         => \$outDir,
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
	print STDERR "OUTPUT DIR: $outDir\n";	
}


my ($counter, @temp, $gene, $score, $type, $threshold, %reportedGenes);

# SIFT
$counter = 1;
$threshold = 0.05;
%reportedGenes = ();
$file_out = "$outDir/SIFT.result";
open(IN, "cut -f 7,13 $file_in | tail -n+2 | sort -k2,2n |");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[1];
	next if ($score eq "." || exists $reportedGenes{$gene});
	last if($score > $threshold);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
	$reportedGenes{$gene} = "";
	$counter++;
}
close(OUT);
close(IN);

# Polyphen2_HVAR_score
$counter = 1;
$threshold = 0.909; # stringent mode
# $threshold = 0.446; # relaxed mode
%reportedGenes = ();
$file_out = "$outDir/Polyphen2_HVAR.result";
open(IN, "cut -f 7,17 $file_in | tail -n+2 | sort -k2,2nr |");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[1];
	next if ($score eq "." || exists $reportedGenes{$gene});
	last if($score < $threshold);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
	$reportedGenes{$gene} = "";
	$counter++;
}
close(OUT);
close(IN);

# MutationAssessor
$counter = 1;
%reportedGenes = ();
$file_out = "$outDir/MutationAssessor.result";
open(IN, "cut -f 7,23,24 $file_in | tail -n+2 | sort -k2,2nr |");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[1];
	$type = $temp[2];
	next if ($score eq "." || exists $reportedGenes{$gene});
	next unless ($type eq "H"	# stringent mode
				 #|| $type eq "M" # relaxed mode
	);
	print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
	$reportedGenes{$gene} = "";
	$counter++;
}
close(OUT);
close(IN);
