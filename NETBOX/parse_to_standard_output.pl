	#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($file_in, $file_out, $mutation_frequency_file, $out_dir, $flag_debug, $flag_help);

my $help_message = "
This script parses NetBox's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to NetBox modules results file *
	--mutation = path to mutation frequency file *
	--outDir = path to output directory *
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
	"mutation=s"	=> \$mutation_frequency_file,
	"outDir=s"      => \$out_dir,
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
	print STDERR "MUTATION: $mutation_frequency_file\n";
	print STDERR "OUTPUT: $out_dir\n";	
}


my %gene_freq = ();
open(FILE, $mutation_frequency_file);
while(<FILE>){
    @line = split(/\t/, $_);
    $gene = $line[0];
    $freq = $line[1];
    $gene_freq{$gene} = $freq;
}
close(FILE);

$file_out = "$out_dir/NetBox.result";
open(OUT, ">$file_out.temp");
open(FILE, $file_in);
<FILE>;#To skip the header
while(<FILE>){
    @line = split(/\s+/, $_);
    $gene = $line[0];
    print OUT $gene."\t".$gene_freq{$gene}."\n";
}
close(FILE);
close(OUT);

open(OUT, ">$file_out");
print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";
open(FILE, "sort -k2,2 -gr $file_out.temp |");
my $rank = 1;
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene = $line[0];
    $score = $line[1];

    print OUT $gene."\t"."ALL"."\t".$rank."\t".$score."\t"."-"."\n";
    $rank++;
}
close(FILE);
close(OUT);

system("rm -f $file_out.temp");