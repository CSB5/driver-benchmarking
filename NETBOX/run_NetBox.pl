#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs NetBox.

Usage:
	run_NetBox.pl [OPTIONS]

Options:
	--config = path to config file *
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

# Preparing output folder
$command = "mkdir $config{'default.outDir'}";
print STDERR "$command\n" if $flag_debug;
system($command) unless (-d $config{'default.outDir'});

open(FILE, $config{'default.gene'});
my %gene_set = ();
while(<FILE>){
    $gene = $_;chop $gene;
    $gene_set{$gene} = 1;
}
close(FILE);

#Get the gene list
open(OUT, ">$config{'default.outDir'}/gene_list.txt");
open(FILE, "sort -k2,2 -gr $config{'default.mutationFrequency'} |");
$nb_selected = 0;
while(<FILE>){
    @line = split(/\t/, $_);
    $gene = $line[0];
    if(exists $gene_set{$gene}){
	$nb_selected++;
	print OUT $gene."\n";
    }
    last if($nb_selected == $config{'default.maxMutation'});
}
close(OUT);
close(FILE);

#Construt the config file
open(OUT, ">$config{'default.outDir'}/netbox1.props");
print OUT "gene_file=$config{'default.outDir'}/gene_list.txt
title=$config{'default.maxMutation'} most freqently mutated genes
shortest_path_threshold=2
p_value_threshold=0.05
num_global_trials=0
num_local_trials=0\n";
close(OUT);

#Run NetBox
system("cd $config{'default.outDir'}; netAnalyze.py netbox1.props");

