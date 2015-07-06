#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($in_dir, $out_dir, $flag_debug, $flag_help);

my $help_message = "
This script parses ActiveDriver's output to a standard output.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--inDir = path to ActiveDriver results directory *
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
	"inDir=s"      	=> \$in_dir,
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
	print STDERR "INPUT: $in_dir\n";
	print STDERR "OUTPUT: $out_dir\n";
}

system("cp $in_dir/*.result $out_dir");
