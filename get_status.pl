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
	1.2

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
my $date = localtime;
print STATUS "Date: $date\n";

# CHASM
print STATUS "CHASM: ";
if(-s "$config{'general.analysisDir'}/CHASM/LATEST/Variant_Additional_Details.Result.tsv"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# DawnRank
print STATUS "DawnRank: ";
if(-s "$config{'general.analysisDir'}/DAWNRANK/LATEST/driver_list.dat"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# DriverNet
print STATUS "DriverNet: ";
if(-s "$config{'general.analysisDir'}/DRIVERNET/LATEST/res_driver_net.dat"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# HotNet2
print STATUS "HotNet2: ";
if(-s "$config{'general.analysisDir'}/HOTNET2/LATEST/HotNet2.result"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# LJB
print STATUS "LJB: ";
if(-s "$config{'LJB.annotation'}"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# MutSigCV
print STATUS "MutSigCV: ";
if(-s "$config{'general.analysisDir'}/MUTSIGCV/LATEST/$config{'general.disease'}.sig_genes.txt"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# NetBox
print STATUS "NetBox: ";
if(-s "$config{'general.analysisDir'}/NETBOX/LATEST/modules.txt"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# OncodriveCIS
print STATUS "OncodriveCIS: ";
if(-s "$config{'general.analysisDir'}/ONCODRIVECIS/LATEST/OncoCNA.combined.txt"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# OncodriveCLUST
print STATUS "OncodriveCLUST: ";
if(-s "$config{'general.analysisDir'}/ONCODRIVECLUST/LATEST/oncodriveclust-results.tsv"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# OncodriveFM
print STATUS "OncodriveFM: ";
if(-s "$config{'general.analysisDir'}/ONCODRIVEFM/LATEST/OncodriveFM-genes.tsv"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# oncoIMPACT
print STATUS "oncoIMPACT: ";
if(-s "$config{'general.analysisDir'}/ONCOIMPACT/LATEST/driver_list.txt"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# oncoIMPACT-v1
print STATUS "oncoIMPACT-v1: ";
if(-s "$config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST/driver_list.txt"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}

# S2N
print STATUS "S2N: ";
if(-s "$config{'general.analysisDir'}/S2N/LATEST/S2N.result"){
	print STATUS "OK\n";
} else{
	print STATUS "Failed\n";
}


close(STATUS);