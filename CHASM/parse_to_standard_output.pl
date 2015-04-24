#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($inDir, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses CHASM's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--inDir = path to CHASM results folder *
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
	"inDir=s"      	=> \$inDir,
	"out=s"     	=> \$file_out,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "INPUT DIR: $inDir\n";
	print STDERR "OUTPUT: $file_out\n";	
}

my (%samples, %fdrs, @temp, $geneID, $sampleID, $geneScore, $fdr);

# Read variant level results to get samples
my $threshold = 0.2;	# FDR used in CHASM paper
open(FILE, "$inDir/Variant_Analysis.Result.tsv");
(%samples, %fdrs) = ();
while(<FILE>){
	next if( $_ =~ /^#/ || $_ =~ /^ID/ );
	chomp(@temp = split(/\t/, $_));
	$geneID = $temp[7];
	$sampleID = $temp[6];
	$fdr = $temp[21];
	next if( $fdr eq " " || $fdr > $threshold );
	unless(exists $samples{$geneID}){
		my @list = ();
		my @list1 = ();
		$samples{$geneID} = \@list;
		$fdrs{$geneID} = \@list1;
	}
	push(@{$samples{$geneID}}, $sampleID);
	push(@{$fdrs{$geneID}}, $fdr);
}
close(FILE);

# Generate report
open(FILE, "grep -v \"#\" $inDir/Gene_Level_Analysis.Result.tsv | tail -n+2 | cut -f 1,4 | sort -k2,2g |");
open(OUT, ">$file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
my $rank = 1;
while(<FILE>){
	chomp(@temp = split(/\t/, $_));
	$geneID = $temp[0];
	$geneScore = $temp[1];
	if(exists $samples{$geneID} && $geneScore eq " "){
		print STDERR "ERROR: $geneID has significant pVal but no gene score!\n";
		next;
	}
	next unless(exists $samples{$geneID});
	print OUT $geneID. "\t" . join(";", @{$samples{$geneID}}) . "\t" . $rank . "\t" . $geneScore . "\t" . "fdr:" . join(";", @{$fdrs{$geneID}}) . "\n";
	$rank++;
}
close(OUT);
close(FILE);
