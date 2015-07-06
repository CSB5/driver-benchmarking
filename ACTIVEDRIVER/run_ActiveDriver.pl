#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script runs ActiveDriver.

Usage:
	run_ActiveDriver.pl [OPTIONS]

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

# Create output directory if it doesn't exists
$command = "mkdir $config{'default.outDir'}";
print STDERR "$command\n" if($flag_debug);
system($command) unless (-d $config{'default.outDir'});

# Run ActiveDriver against default database
print "Running ActiveDriver against phospho sites. Please wait...";
$command = "$config{'default.scriptsDir'}/ActiveDriver_phospho.R -a $config{'default.annotation'} -f $config{'default.outDir'}/ActiveDriver_phospho.result";
print STDERR "$command\n" if $flag_debug;
system($command);
print "done.\n";

# Run ActiveDriver against PTM sites
print "Running ActiveDriver against PTM sites. Please wait...";
$command = "$config{'default.scriptsDir'}/ActiveDriver_PTM.R -a $config{'default.annotation'} -f $config{'default.outDir'}/ActiveDriver_PTM.result";
print STDERR "$command\n" if $flag_debug;
system($command);
print "done.\n";
