#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script runs OncodriveFM.

Usage:
	run_OncodriveFM.pl [OPTIONS]

Options:
	--config = path to config file *
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
	"config=s" => \$configFile,
	"debug"    => \$flag_debug,
	"help"     => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

# Read config file
print "Reading config file. Please wait...";
Config::Simple->import_from( $configFile, \%config );
print "done.\n";

my ($annot_file, $out_dir);
$annot_file = $config{'default.annotation'};
$out_dir = $config{'default.outDir'};


# Generate OncodriveFM input
if (index($annot_file, "hg19") != -1) {
	# hg19
	my %score_col = (
	    "SIFT_score", 12,
	    "Polyphen2_HDIV_score", 14,
	    "Polyphen2_HVAR_score", 16,
	    "LRT_score", 18,
	    "MutationTaster_score", 20,
	    "MutationAssessor_score", 22,
	    "FATHMM_score", 24,
	    "RadialSVM_score", 26,
	    "LR_score", 28,
	    "VEST3_score", 30,
    );    
} else{
	# hg18
	my %score_col = (
	    "SIFT_score", 10,
	    "Polyphen2_HDIV_score", 12,
	    "Polyphen2_HVAR_score", 14,
	    "LRT_score", 16,
	    "MutationTaster_score", 18,
	    "MutationAssessor_score", 20,
	    "FATHMM_score", 22,
	    "RadialSVM_score", 24,
	    "LR_score", 26,
	    "VEST3_score", 28,
    );   
}



my %score_used = (
    "SIFT_score", "SIFT",
    "Polyphen2_HVAR_score", "PPH2",
    "MutationAssessor_score", "MA"
    );
my @score_order = keys %score_used;

$header = "SAMPLE\tGENE";
foreach $s (@score_order){
    $header .= "\t".$score_used{$s};
}
open(OUT, ">$out_dir/OncodriveFM.txt");
print OUT $header."\n";

open(FILE, $annot_file);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    
    next if($line[0] eq "Chr");


    $all_define = 1;
    my @res_score = ();
    for($i = 0; $i < @score_order; $i++){
	$s = $score_order[$i];
	$s_val = $line[$score_col{$s}];
	push(@res_score, $s_val);
	if($s_val eq "."){
	    $all_define = 0;
	}
    }
    
    if($all_define){
	$sample_ID = $line[@line-1];
	$gene = $line[6];
	print OUT $sample_ID."\t".$gene."\t".(join("\t", @res_score))."\n";
    }
}

close(FILE);
close(OUT);


# Run OncodriveFM
$command = "export PATH=/mnt/software/unstowable/anaconda/bin/:\$PATH; source activate oncodrivefm; oncodrivefm -e median -m $config{'default.mappingFile'} -j $config{'default.numThreads'} -o $out_dir/ $out_dir/OncodriveFM.txt";
print STDERR "$command\n" if $flag_debug;
system($command);

