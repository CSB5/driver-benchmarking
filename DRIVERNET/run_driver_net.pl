#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $nb_random_test_by_proc, $flag_debug, $flag_help);

my $help_message = "
This script runs DriverNet.

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


$nb_random_test_by_proc = 500 / $config{'default.numProc'};

$driver_net_dir = "$config{'default.outDir'}";
system("mkdir -p $driver_net_dir") unless (-d $driver_net_dir);

	
open (OUT_cmd, ">$driver_net_dir/cmd.txt");
	
#Construct the R file
for($i = 0; $i < $config{'default.numProc'}; $i++){
    $file = "$driver_net_dir/rand_$i.R";
    open(OUT, ">$file");
    print OUT "library(DriverNet)\n";
    print OUT "data(sampleGeneNames)\n";
    print OUT "my_patMut <- read.table(\"$config{'default.mutMatrix'}\", header=T)\n";
    print OUT "load(\"$config{'default.scriptsDir'}/influenceGraph.rda\")\n";

    #Data given in the paper
    if(index($config{'default.expData'}, ".rda") != -1){
	print OUT "load(\"$config{'default.expData'}\")\n";
    }
    #Use an outlier matrix computed from the data
    else{
	print OUT "patExpMatrix <- read.table(\"$config{'default.expData'}\", header=T)\n";
	print OUT "patOutMatrix <- getPatientOutlierMatrix (as.matrix(patExpMatrix), th=2)\n";
    }

    print OUT "randomDriversResult = computeRandomizedResult(patMutMatrix=my_patMut, patOutMatrix=patOutMatrix, influenceGraph=influenceGraph, geneNameList= sampleGeneNames, outputFolder=NULL, printToConsole=FALSE,numberOfRandomTests=$nb_random_test_by_proc, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)\n";
    
    print OUT "save(randomDriversResult, file=\"$driver_net_dir/rand_$i.rda\")\n";
    close(OUT);
	    
    print OUT_cmd "/mnt/software/bin/Rscript-3.1.0 $file &> $file.log\n";
}
close(OUT_cmd);
	
#run the analysis in parallel
$cmd = "cat $driver_net_dir/cmd.txt | xargs -I cmd --max-procs=$config{'default.numProc'} bash -c cmd > /dev/null \n";
system("$cmd");

#Combine the results
$file = "$driver_net_dir/rand_comb.R";
open(OUT, ">$file");
for($i = 0; $i < $config{'default.numProc'}; $i++){
    print OUT "load(\"$driver_net_dir/rand_$i.rda\")\n";
    print OUT "X_$i = randomDriversResult\n";
    print OUT "randomDriversResult_comb = c(X_0, X_1)\n" if($i == 1);
    print OUT "randomDriversResult_comb = c(randomDriversResult_comb, X_$i)\n" if($i > 1);
}
print OUT "save(randomDriversResult_comb, file=\"$driver_net_dir/rand_comb.rda\")\n";
close(OUT);
system("/mnt/software/bin/Rscript-3.1.0 $file\n");

#Get the pvalue
$file = "$driver_net_dir/driver_net.R";
open(OUT, ">$file");
print OUT "library(DriverNet)\n";
print OUT "data(sampleGeneNames)\n";
print OUT "my_patMut <- read.table(\"$config{'default.mutMatrix'}\", header=T)\n";
print OUT "load(\"$config{'default.scriptsDir'}/influenceGraph.rda\")\n";
#
#Data given in the paper
if(index($config{'default.expData'}, ".rda") != -1){
    print OUT "load(\"$config{'default.expData'}\")\n";
}
#Use an outlier matrix computed from the data
else{
    print OUT "patExpMatrix <- read.table(\"$config{'default.expData'}\", header=T)\n";
    print OUT "patOutMatrix <- getPatientOutlierMatrix(as.matrix(patExpMatrix), th=2)\n";
    
}
#
#
print OUT "driversList = computeDrivers(my_patMut,patOutMatrix,influenceGraph,outputFolder=NULL, printToConsole=FALSE)\n";
print OUT "load(\"$driver_net_dir/rand_comb.rda\")\n";
print OUT "res = resultSummary(driversList, randomDriversResult_comb, my_patMut, influenceGraph, outputFolder=NULL, printToConsole=FALSE)\n";
print OUT "write.table(res[,1:4], file=\"$driver_net_dir/res_driver_net.dat\", sep=\"\t\", row.names=T, col.names=T, quote=FALSE)\n";
close(OUT);
system("/mnt/software/bin/Rscript-3.1.0 $file\n");
    
