#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ( $samplesList, $TCGAmaf, $outDir, $flag_debug, $flag_help );

my $help_message = "
This script filters the TCGA Somatic Mutations MAF file against a list of samples.

Usage:
	filter_maf.pl [OPTIONS]

Options:
	--samples = path to file containing list of samples *
	--maf = path to TCGA somatic mutations MAF file *
	--outDir = path to output directory
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
	"samples=s" => \$samplesList,
	"maf=s"     => \$TCGAmaf,
	"outDir=s"  => \$outDir,
	"debug"     => \$flag_debug,
	"help"      => \$flag_help
  )
  or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[DEBUG] Samples:$samplesList\n";
	print STDERR "[DEBUG] TCGA MAF:$TCGAmaf\n";
	print STDERR "[DEBUG] Output directory:$outDir\n";
}

my ( %samples, $sampleID, @temp );

# Reading list of samples
open( FILE, $samplesList );
while (<FILE>) {
	$samples{$_} = "";
}
close(FILE);

# Generating filtered MAF
open( OUT,  "> $outDir/TCGA_somatic_mutations.filtered.maf" );
print STDERR "[DEBUG] Output file:$outDir/TCGA_somatic_mutations.filtered.maf\n" if $flag_debug; 
open( FILE, $TCGAmaf );
print OUT <FILE>;    # print header
while (<FILE>) {
	chomp( @temp = split( /\t/, $_ ) );
	$sampleID = $temp[15];
	@temp     = split( /-/, $sampleID );
	$sampleID = "$temp[0]-$temp[1]-$temp[2]-$temp[3]";
	print STDERR "Checking sample:$sampleID\n" if $flag_debug;
	print OUT $_ if ( exists $samples{$sampleID} );
}
close(OUT);
close(FILE);
