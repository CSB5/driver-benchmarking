#!/usr/bin/perl

use warnings;
use Getopt::Long;
use Config::Simple;
use POSIX 'strftime';

my ($configFile, $flag_help, %config);

my $help_message = "
This script runs produces the status information of various softwares for benchmarking.

Usage:
	get_status.pl [OPTIONS]

Options:
	--config = path to config file *
	--help : prints this message 
	
* indicates required parameters	


Version:
	$version

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}

GetOptions(
	"config=s" => \$configFile,
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


open( STATUS, "> $config{'general.analysisDir'}/LOGS/status.log" );
my $date = strftime '%Y%m%d', localtime;
print STATUS "Date: $date\n";

# oncoIMPACT
print STATUS "oncoIMPACT: ";
if(-s "$config{'general.analysisDir'}/ONCOIMPACT/LATEST/ANALYSIS/GENE_LIST/ALTERATION.dat"){
	print "OK\n";
} else{
	print "Failed\n";
}

# DriverNet
print STATUS "DriverNet: ";
if(-s "$config{'general.analysisDir'}/DRIVERNET/LATEST/res_driver_net.dat"){
	print "OK\n";
} else{
	print "Failed\n";
}

#MutSigCV
print STATUS "MutSigCV: ";
if(-s "$config{'general.analysisDir'}/MUTSIGCV/LATEST/$config{'general.disease'}.sig_genes.txt"){
	print "OK\n";
} else{
	print "Failed\n";
}

# DawnRank
print STATUS "DawnRank: ";
if(-s "$config{'general.analysisDir'}/DAWNRANK/LATEST/driver_list.dat"){
	print "OK\n";
} else{
	print "Failed\n";
}

# OncodriveFM
print STATUS "OncodriveFM: ";
if(-s "$config{'general.analysisDir'}/ONCODRIVEFM/LATEST/OncodriveFM-genes.tsv"){
	print "OK\n";
} else{
	print "Failed\n";
}

# LJB
print STATUS "LJB: ";
if(-s "$config{'LJB.annotation'}"){
	print "OK\n";
} else{
	print "Failed\n";
}

# OncodriveCLUST
print STATUS "OncodriveCLUST: ";
if(-s "$config{'general.analysisDir'}/ONCODRIVECLUST/LATEST/oncodriveclust-results.tsv"){
	print "OK\n";
} else{
	print "Failed\n";
}

# NetBox
print STATUS "NetBox: ";
if(-s "$config{'general.analysisDir'}/NETBOX/LATEST/modules.txt"){
	print "OK\n";
} else{
	print "Failed\n";
}

close(STATUS);