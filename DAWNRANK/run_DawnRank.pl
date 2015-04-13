#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs DawnRank.

Usage:
	run_DawnRank.pl [OPTIONS]

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
system($command) unless (-e $config{'default.outDir'});


# constructing network matrix required by DawnRank
print "Constructing network matrix. Please wait...";
$command = "$config{'default.scriptsDir'}/construct_adj_matrix.pl --adj $config{'default.adj'} --out $config{'default.outDir'}/network.dat";
$command = $command . " --debug" if $flag_debug;
print STDERR "[command] $command\n" if $flag_debug;
system($command);
print "done.\n";


# filter and sort expression and mutation matrices
print "Filtering and sorting expression and mutation matrices. Please wait...";
$command = "$config{'default.scriptsDir'}/filter_and_update_matrices.pl --exp $config{'default.exp'} --mut $config{'default.mut'} --adj $config{'default.adj'} --out $config{'default.outDir'}";
$command = $command . " --debug" if $flag_debug;
print STDERR "[command] $command\n" if $flag_debug;
system($command);
print "done.\n";


# Running DawnRank
print "Running DawnRank. Please wait...";
$command = "$config{'default.scriptsDir'}/DawnRank.R $config{'default.outDir'}/expression.dat $config{'default.outDir'}/mutation.dat $config{'default.outDir'}/network.dat $config{'default.outDir'}";
print STDERR "[command] $command\n" if $flag_debug;
system($command);
print "done.\n";
