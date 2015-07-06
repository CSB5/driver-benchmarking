#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);
my $qsub = "qsub -terse -m ae -M \$USER_PRINCIPAL_NAME -cwd -V";

my $help_message = "
This script runs HotNet2.

Usage:
	run_HotNet2.pl [OPTIONS]

Options:
	--config = path to config file *
	--debug: prints trace to STDERR
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

my ($geneMutationFile, $irefRunFolder, $irefConfig, $irefResults, $jid);
my $runtime = "$config{'default.runtime'}:0:0";

# Preparing output folder
$command = "mkdir $config{'default.outDir'}";
system($command) unless (-e $config{'default.outDir'});


# Preprocess gene mutation file
$geneMutationFile = "$config{'default.outDir'}/gene_mutation_frequency.dat";
$command = "tail -n+2 $config{'default.freqFile'} | cut -f 1,2 > $geneMutationFile";
system($command);
print STDERR "$command\n" if ($flag_debug);


# Generate scripts
$irefRunFolder = "$config{'default.outDir'}/IREF";
$irefConfig = "$config{'default.outDir'}/iref.cfg";

## Generating irefConfig
open(OUT, ">$irefConfig");
print OUT "--runname iref\n";
print OUT "--infmat_file $config{'default.irefMaf'}\n";
print OUT "--infmat_index_file $config{'default.irefIndex'}\n";
print OUT "--permuted_networks_path $config{'default.irefNetworks'}\n";
print OUT "--heat_file $geneMutationFile\n";
print OUT "--output_directory $irefRunFolder\n";
print OUT "--delta_permutations $config{'default.deltaPerm'}\n";
print OUT "--significance_permutations $config{'default.sigPerm'}\n";
print OUT "--num_cores $config{'default.numThreads'}\n";
close OUT;


# iref run
## run HotNet2
$command = "$qsub -l mem_free=$config{'default.mem'}G,h_rt=$runtime -pe OpenMP $config{'default.numThreads'} -N runHotNet2_iref $config{'default.installationDir'}/runHotNet2.py \@$irefConfig";
print STDERR "$command\n" if ($flag_debug);
$jid = `$command`;
chomp($jid);

## chooseDelta
$command = "$qsub -l mem_free=$config{'default.mem'}G,h_rt=$runtime -pe OpenMP 1 -hold_jid $jid -N chooseDelta_iref $config{'default.scriptsDir'}/choseDelta-v2.py $irefRunFolder $geneMutationFile $config{'default.resultsDir'} $config{'default.sigThreshold'}";
print STDERR "$command\n" if ($flag_debug);
system($command);
