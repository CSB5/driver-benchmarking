#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config);

my $help_message = "
This script runs MutSigCV.

Usage:
	run_driver_net.pl [OPTIONS]

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

my $command;

$command = "cd $config{'default.outDir'};";	# set working directory to outDir
$command = $command . "run_MutSigCV.sh $config{'default.matlab'} $config{'default.maf'} $config{'default.coverage'} $config{'default.covariate'} $config{'default.prefix'} $config{'default.dict'} $config{'default.chr'}";	# run MutSigCV

print STDERR "[MutSigCV] Command:$command\n" if $flag_debug;
#system("ssh -t aries \"$command\"");
system($command);

