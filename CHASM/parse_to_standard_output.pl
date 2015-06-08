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
	1.1

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

my (%samples, %fdrs, %bestScore, @temp, $geneID, $sampleID, $geneScore, $fdr);

# Read variant level results to get samples
my $threshold = 0.2;	# FDR used in CHASM paper
open(FILE, "$inDir/Variant_Additional_Details.Result.tsv");
(%samples, %fdrs, %bestScore) = ();
while(<FILE>){
	next if( $_ =~ /^#/ || $_ =~ /^Input line/ || $_ eq "\n");
	chomp(@temp = split(/\t/, $_));
	$geneID = $temp[8];
	$sampleID = $temp[7];
	$geneScore = $temp[15];
	$fdr = $temp[17];	
	next if( $fdr eq "" || $fdr > $threshold );
	unless(exists $samples{$geneID}){
		my @list = ();
		$samples{$geneID} = \@list;
		$fdrs{$geneID} = $fdr;
		$bestScore{$geneID} = $geneScore;
	} 
	push(@{$samples{$geneID}}, $sampleID);
	if($geneScore > $bestScore{$geneID}){
		$bestScore{$geneID} = $geneScore; 
		$fdrs{$geneID} = $fdr;		
	}
}
close(FILE);

# Generate report
open(OUT, ">$file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
my $rank = 1;
foreach $geneID (sort { $bestScore{$b} <=> $bestScore{$a} or $a cmp $b } keys %bestScore) {
	print OUT $geneID. "\t" . join(";", @{$samples{$geneID}}) . "\t" . $rank . "\t" . $bestScore{$geneID} . "\t" . "fdr:" . $fdrs{$geneID} . "\n";
	$rank++;
}
close(OUT);
close(FILE);
