#!/usr/bin/perl

use warnings;
use Getopt::Long;
use File::Basename;

my ($file_out, $dir, $flag_help, $flag_debug);

my $help_message = "
This script parses oncoIMPACT's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to oncoIMPACT results folder *
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
	"in=s"			=> \$dir,
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
	print STDERR "DIRECTORY: $dir\n";
	print STDERR "OUTPUT: $file_out\n";	
}


my ($counter, @temp, $gene, $score, %ht, @samples, $sampleID);

# Read OncoIMPACT output directory to get list of samples
opendir( DIR, "$dir/stringent/samples" );
# Populate hash with mutation information
while (my $file = readdir(DIR)) {
	next unless (-f "$dir/stringent/samples/$file");
	next unless ($file =~ m/\.tsv$/);
	
	$sampleID = basename($file, ".tsv");
	open(FILE, "$dir/stringent/samples/$file") or die ("File not found: $dir/stringent/samples/$file");
	<FILE>; #skip header
	while(<FILE>){		
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		unless(exists $ht{$gene}){
			my @sampleList = ();
			$ht{$gene} = \@sampleList;
		}
		push( @{$ht{$gene}}, $sampleID);
	}
	close(FILE);
}
close(DIR);

# Generate results
my $score_threshold = 0;
$counter = 1;
open(IN, "$dir/driver_list.txt");
open(OUT, "> $file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
<IN>;	# skip header
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	$score = $temp[7];
	last if($score <= $score_threshold);
	print OUT $gene . "\t" . join(";", @{$ht{$gene}}) . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
	$counter++;
}
close(OUT);
close(IN);
