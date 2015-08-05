#!/usr/bin/perl

#use warnings;
#no warnings 'experimental';
use Config::Simple;
use Getopt::Long;
use POSIX 'strftime';
use 5.010;

my $version = "v3.9.1";
my $date = strftime '%Y%m%d', localtime;
my $runID = "${date}_${version}";

my ( $configFile, $flag_debug, $flag_help, $flag_update, $flag_simulate, %config, @queue, $command, $lastID, $software, $outDir, $mutSigID, $oncoIMPACT_ID);
my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V";

my $help_message = "
This script runs various softwares for benchmarking.

Usage:
	run_benchmark_pipeline.pl [OPTIONS]

Options:
	--config = path to config file *
	--debug: prints trace to STDERR
	--sim: simulation mode - only prints commands to be executed to TRACE log
	--update: update mode - only regenerates consolidated results
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
	"debug"    => \$flag_debug,
	"update"   => \$flag_update,
	"sim"	   => \$flag_simulate,
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


# Update only mode
if($flag_update){

	print "Updating latest consolidated results. Please wait.\n";
	my $resultsDir = "$config{'general.analysisDir'}/CONSOLIDATED_RESULTS/LATEST";
	my $logsDir = "$config{'general.analysisDir'}/LOGS/$runID";
	my $numSamples;
	$numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);

	# ActiveDriver
	print "ActiveDriver: ";
	$command = "$config{'general.scriptsDir'}/ACTIVEDRIVER/parse_to_standard_output.pl --inDir $config{'general.analysisDir'}/ACTIVEDRIVER/LATEST --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/ACTIVEDRIVER/LATEST/ActiveDriver_phospho.result" && -s "$config{'general.analysisDir'}/ACTIVEDRIVER/LATEST/ActiveDriver_PTM.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# CHASM
	print "CHASM: ";
	$command = "$config{'general.scriptsDir'}/CHASM/parse_to_standard_output.pl --inDir $config{'general.analysisDir'}/CHASM/LATEST --samples $numSamples --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/CHASM/LATEST/Variant_Additional_Details.Result.tsv"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# DawnRank
	print "DawnRank: ";
	$command = "$config{'general.scriptsDir'}/DAWNRANK/parse_to_standard_output.pl --in $config{'general.analysisDir'}/DAWNRANK/LATEST/driver_list.dat --out $resultsDir/DawnRank.result";
	if(-s "$config{'general.analysisDir'}/DAWNRANK/LATEST/driver_list.dat"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# DriverNet
	print "DriverNet: ";
	$command = "$config{'general.scriptsDir'}/DRIVERNET/parse_to_standard_output.pl --in $config{'general.analysisDir'}/DRIVERNET/LATEST/res_driver_net.dat --out $resultsDir/DriverNet.result";
	if(-s "$config{'general.analysisDir'}/DRIVERNET/LATEST/res_driver_net.dat"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# FATHMM
	print "FATHMM: ";
	$command = "$config{'general.scriptsDir'}/FATHMM/parse_to_standard_output.pl --in $config{'general.analysisDir'}/FATHMM/LATEST/fathmm.result --samples $numSamples --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/FATHMM/LATEST/fathmm.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# HotNet2
	print "HotNet2: ";
	$command = "$config{'general.scriptsDir'}/HOTNET2/parse_to_standard_output.pl --in $config{'general.analysisDir'}/HOTNET2/LATEST/HotNet2.result --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/HOTNET2/LATEST/HotNet2.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# HotNet2 alternative
	print "HotNet2A: ";
	$command = "$config{'general.scriptsDir'}/HOTNET2A/parse_to_standard_output.pl --in $config{'general.analysisDir'}/HOTNET2/LATEST/HotNet2.result --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/HOTNET2A/LATEST/HotNet2.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# LJB
	print "LJB: ";
	$command = "$config{'general.scriptsDir'}/LJB/parse_to_standard_output.pl --in $config{'LJB.annotation'} --samples $numSamples --outDir $resultsDir";
	if(-s "$config{'LJB.annotation'}"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# MutSigCV
	print "MutSigCV: ";
	$command = "$config{'general.scriptsDir'}/MUTSIGCV/parse_to_standard_output.pl --in $config{'general.analysisDir'}/MUTSIGCV/LATEST/$config{'general.disease'}.sig_genes.txt --out $resultsDir/MutSigCV.result";
	if(-s "$config{'general.analysisDir'}/MUTSIGCV/LATEST/$config{'general.disease'}.sig_genes.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# NetBox
	print "NetBox: ";
	$command = "$config{'general.scriptsDir'}/NETBOX/parse_to_standard_output.pl --in $config{'general.analysisDir'}/NETBOX/LATEST/modules.txt --mutation $config{'NetBox.mutationFrequency'} --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/NETBOX/LATEST/modules.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncodriveCIS
	print "OncodriveCIS: ";
	$command = "$config{'general.scriptsDir'}/ONCODRIVECIS/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCODRIVECIS/LATEST/OncoCNA.combined.txt --out $resultsDir/OncodriveCIS.result";
	if(-s "$config{'general.analysisDir'}/ONCODRIVECIS/LATEST/OncoCNA.combined.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncodriveCLUST
	print "OncodriveCLUST: ";
	$command = "$config{'general.scriptsDir'}/ONCODRIVECLUST/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCODRIVECLUST/LATEST/oncodriveclust-results.tsv --out $resultsDir/OncodriveCLUST.result";
	if(-s "$config{'general.analysisDir'}/ONCODRIVECLUST/LATEST/oncodriveclust-results.tsv"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncodriveFM
	print "OncodriveFM: ";
	$command = "$config{'general.scriptsDir'}/ONCODRIVEFM/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCODRIVEFM/LATEST/OncodriveFM-genes.tsv --out $resultsDir/OncodriveFM.result";
	if(-s "$config{'general.analysisDir'}/ONCODRIVEFM/LATEST/OncodriveFM-genes.tsv"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncoIMPACT
	print "OncoIMPACT: ";
	$command = "$config{'general.scriptsDir'}/ONCOIMPACT/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCOIMPACT/LATEST/ --out $resultsDir/oncoIMPACT.result";
	if(-s "$config{'general.analysisDir'}/ONCOIMPACT/LATEST/driver_list.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncoIMPACT
	print "oncoIMPACT-discovery: ";
	$command = "$config{'general.scriptsDir'}/ONCOIMPACT/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCOIMPACT-DISCOVERY/LATEST/ --out $resultsDir/oncoIMPACT-discovery.result";
	if(-s "$config{'general.analysisDir'}/ONCOIMPACT-DISCOVERY/LATEST/driver_list.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# OncoIMPACT-v1
	print "OncoIMPACT-v1: ";
	$command = "$config{'general.scriptsDir'}/ONCOIMPACT-V1/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST/ --out $resultsDir/oncoIMPACT-v1.result";
	if(-s "$config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST/driver_list.txt"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# S2N
	print "S2N: ";
	$command = "$config{'general.scriptsDir'}/S2N/parse_to_standard_output.pl --in $config{'general.analysisDir'}/S2N/LATEST/S2N.result --out $resultsDir/S2N.result";
	if(-s "$config{'general.analysisDir'}/S2N/LATEST/S2N.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	# transFIC
	print "transFIC: ";
	$command = "$config{'general.scriptsDir'}/TRANSFIC/parse_to_standard_output.pl --in $config{'general.analysisDir'}/TRANSFIC/LATEST/transFIC.result --samples $numSamples --outDir $resultsDir";
	if(-s "$config{'general.analysisDir'}/TRANSFIC/LATEST/transFIC.result"){
		system($command);
		print "Done\n";
	} else{
		print "Incomplete run!\n";
	}

	exit 0;
}


# Preparing output folders
print "Preparing output folders. Please wait...";
my $analysisDir;
my $resultsDir = "$config{'general.analysisDir'}/CONSOLIDATED_RESULTS/$runID";
unless (-e $resultsDir){
	system("mkdir -p $resultsDir");
	system("cp -p $config{'general.analysisDir'}/CONSOLIDATED_RESULTS/LATEST/* $resultsDir") if(-e "$config{'general.analysisDir'}/CONSOLIDATED_RESULTS/LATEST/");
	system("ln -sfn $resultsDir $config{'general.analysisDir'}/CONSOLIDATED_RESULTS/LATEST");
	system("chmod g+w $resultsDir");
}
my $logsDir = "$config{'general.analysisDir'}/LOGS/$runID";
system("mkdir -p $logsDir") unless (-e $logsDir);
print "done.\n";

# Initializing variables
print "Initializing variables. Please wait...";
open( TRACE, ">> $logsDir/trace.log" );
@queue = ();
my $runtime = "$config{'cluster.runtime'}:0:0";
print "done.\n";



### Softwares/Packages ###
# oncoIMPACT
if($config{'general.oncoIMPACT'}){
	print "Running oncoIMPACT. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCOIMPACT/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCOIMPACT/LATEST");
	}

	# Generate config file
	generateConfig("oncoIMPACT");

	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_oncoIMPACT -e $logsDir/oncoIMPACT.error.log -o $logsDir/oncoIMPACT.run.log $config{'general.scriptsDir'}/ONCOIMPACT/oncoIMPACT.pl $analysisDir/oncoIMPACT_$runID.cfg";
	$command = $command . " 1" if ($flag_debug);
	submit($command);

	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$oncoIMPACT_ID = $lastID;
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT_parseOutput -e $logsDir/oncoIMPACT_parseOutput.error.log -o $logsDir/oncoIMPACT_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCOIMPACT/parse_to_standard_output.pl --in $analysisDir --out $resultsDir/oncoIMPACT.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# DriverNet
if($config{'general.DriverNet'}){
	print "Running DriverNet. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/DRIVERNET/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/DRIVERNET/LATEST");
	}

	# Generate config file
	generateConfig("DriverNet");

	# Run DriverNet
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_DriverNet -e $logsDir/DriverNet.error.log -o $logsDir/DriverNet.run.log $config{'general.scriptsDir'}/DRIVERNET/run_driver_net.pl --config $analysisDir/DriverNet_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse DriverNet output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_DriverNet_parseOutput -e $logsDir/DriverNet_parseOutput.error.log -o $logsDir/DriverNet_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/DRIVERNET/parse_to_standard_output.pl --in $analysisDir/res_driver_net.dat --out $resultsDir/DriverNet.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# MutSigCV
if($config{'general.MutSigCV'}){
	print "Running MutSigCV. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/MUTSIGCV/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/MUTSIGCV/LATEST");
	}

	# Generate config file
	generateConfig("MutSigCV");

	# Run MutSigCV
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV -e $logsDir/MutSigCV.error.log -o $logsDir/MutSigCV.run.log $config{'general.scriptsDir'}/MUTSIGCV/run_MutSigCV.pl --config $analysisDir/MutSigCV_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse MutSigCV output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$mutSigID = $lastID; # ID to control HotNet2A run
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV_parseOutput -e $logsDir/MutSigCV_parseOutput.error.log -o $logsDir/MutSigCV_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/MUTSIGCV/parse_to_standard_output.pl --in $analysisDir/$config{'general.disease'}.sig_genes.txt --out $resultsDir/MutSigCV.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# OncodriveFM
if($config{'general.OncodriveFM'}){
	print "Running OncodriveFM. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCODRIVEFM/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCODRIVEFM/LATEST");
	}

	# Generate config file
	generateConfig("OncodriveFM");

	# Run OncodriveFM
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_OncodriveFM -e $logsDir/OncodriveFM.error.log -o $logsDir/OncodriveFM.run.log $config{'general.scriptsDir'}/ONCODRIVEFM/run_OncodriveFM.pl --config $analysisDir/OncodriveFM_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse OncodriveFM output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_OncodriveFM_parseOutput -e $logsDir/OncodriveFM_parseOutput.error.log -o $logsDir/OncodriveFM_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCODRIVEFM/parse_to_standard_output.pl --in $analysisDir/OncodriveFM-genes.tsv --out $resultsDir/OncodriveFM.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# LJB
if($config{'general.LJB'}){
	print "Running LJB. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/LJB/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/LJB/LATEST");
	}

	my $numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_LJB_parseOutput -e $logsDir/LJB_parseOutput.error.log -o $logsDir/LJB_parseOutput.run.log $config{'general.scriptsDir'}/LJB/parse_to_standard_output.pl --in $config{'LJB.annotation'} --samples $numSamples --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# OncodriveCLUST
if($config{'general.OncodriveCLUST'}){
	print "Running OncodriveCLUST. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCODRIVECLUST/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCODRIVECLUST/LATEST");
	}

	# Generate config file
	generateConfig("OncodriveCLUST");

	# Run OncodriveCLUST
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_OncodriveCLUST -e $logsDir/OncodriveCLUST.error.log -o $logsDir/OncodriveCLUST.run.log $config{'general.scriptsDir'}/ONCODRIVECLUST/run_OncodriveCLUST.pl --config $analysisDir/OncodriveCLUST_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse OncodriveCLUST output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_OncodriveCLUST_parseOutput -e $logsDir/OncodriveCLUST_parseOutput.error.log -o $logsDir/OncodriveCLUST_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCODRIVECLUST/parse_to_standard_output.pl --in $analysisDir/oncodriveclust-results.tsv --out $resultsDir/OncodriveCLUST.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# DawnRank
if($config{'general.DawnRank'}){
	print "Running DawnRank. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/DAWNRANK/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/DAWNRANK/LATEST");

	# Generate config file
	generateConfig("DawnRank");

	# Run DawnRank
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=700:0:0 -pe OpenMP 1 -N $config{'general.disease'}_DawnRank -e $logsDir/DawnRank.error.log -o $logsDir/DawnRank.run.log $config{'general.scriptsDir'}/DAWNRANK/run_DawnRank.pl --config $analysisDir/DawnRank_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse DawnRank output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_DawnRank_parseOutput -e $logsDir/DawnRank_parseOutput.error.log -o $logsDir/DawnRank_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/DAWNRANK/parse_to_standard_output.pl --in $analysisDir/driver_list.dat --out $resultsDir/DawnRank.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# NetBox
if($config{'general.NetBox'}){
	print "Running NetBox. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/NETBOX/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/NETBOX/LATEST");

	# Generate config file
	generateConfig("NetBox");

	# Run NetBox
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_NetBox -e $logsDir/NetBox.error.log -o $logsDir/NetBox.run.log $config{'general.scriptsDir'}/NETBOX/run_NetBox.pl --config $analysisDir/NetBox_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse NetBox output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_NetBox_parseOutput -e $logsDir/NetBox_parseOutput.error.log -o $logsDir/NetBox_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/NETBOX/parse_to_standard_output.pl --in $analysisDir/modules.txt --mutation $config{'NetBox.mutationFrequency'} --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# CHASM
if($config{'general.CHASM'}){
	print "Running CHASM. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/CHASM/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/CHASM/LATEST");

	# Generate config file
	generateConfig("CHASM");

	# Run CHASM
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_CHASM -e $logsDir/CHASM.error.log -o $logsDir/CHASM.run.log $config{'general.scriptsDir'}/CHASM/run_CHASM.pl --config $analysisDir/CHASM_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse CHASM output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	my $numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_CHASM_parseOutput -e $logsDir/CHASM_parseOutput.error.log -o $logsDir/CHASM_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/CHASM/parse_to_standard_output.pl --inDir $analysisDir --samples $numSamples --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# oncoIMPACT-v1
if($config{'general.oncoIMPACT-v1'}){
	print "Running oncoIMPACT-v1. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCOIMPACT-V1/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST");
	}

	# Generate config file
	generateConfig("oncoIMPACT-v1");

	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_oncoIMPACT-v1 -e $logsDir/oncoIMPACT-v1.error.log -o $logsDir/oncoIMPACT-v1.run.log -b y $config{'general.scriptsDir'}/ONCOIMPACT-V1/oncoIMPACT.exe --database $analysisDir/oncoIMPACT-v1_$runID.cfg";
	submit($command);

	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT-v1_parseOutput -e $logsDir/oncoIMPACT-v1_parseOutput.error.log -o $logsDir/oncoIMPACT-v1_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCOIMPACT-V1/parse_to_standard_output.pl --in $analysisDir --out $resultsDir/oncoIMPACT-v1.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# HotNet2
if($config{'general.HotNet2'}){
	print "Running HotNet2. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/HOTNET2/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/HOTNET2/LATEST");

	# Generate config file
	generateConfig("HotNet2");

	# Run HotNet2
	$command = "$qsub -l mem_free=200G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_HotNet2 -e $logsDir/HotNet2.error.log -o $logsDir/HotNet2.run.log $config{'general.scriptsDir'}/HOTNET2/run_HotNet2.pl --config $analysisDir/HotNet2_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse HotNet2 output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_HotNet2_parseOutput -e $logsDir/HotNet2_parseOutput.error.log -o $logsDir/HotNet2_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/HOTNET2/parse_to_standard_output.pl --in $analysisDir/HotNet2.result --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# S2N
if($config{'general.S2N'}){
	print "Running S2N. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/S2N/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/S2N/LATEST");

	# Generate config file
	generateConfig("S2N");

	# Run S2N
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_S2N -e $logsDir/S2N.error.log -o $logsDir/S2N.run.log $config{'general.scriptsDir'}/S2N/run_S2N.pl --config $analysisDir/S2N_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse S2N output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_S2N_parseOutput -e $logsDir/S2N_parseOutput.error.log -o $logsDir/S2N_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/S2N/parse_to_standard_output.pl --in $analysisDir/S2N.result --out $resultsDir/S2N.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# OncodriveCIS
if($config{'general.OncodriveCIS'}){
	print "Running OncodriveCIS. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCODRIVECIS/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCODRIVECIS/LATEST");

	# Generate config file
	generateConfig("OncodriveCIS");

	# Run OncodriveCIS
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_OncodriveCIS -e $logsDir/OncodriveCIS.error.log -o $logsDir/OncodriveCIS.run.log $config{'general.scriptsDir'}/ONCODRIVECIS/run_OncodriveCIS.pl --config $analysisDir/OncodriveCIS_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse OncodriveCIS output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_OncodriveCIS_parseOutput -e $logsDir/OncodriveCIS_parseOutput.error.log -o $logsDir/OncodriveCIS_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCODRIVECIS/parse_to_standard_output.pl --in $analysisDir/OncoCNA.combined.txt --out $resultsDir/OncodriveCIS.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# ActiveDriver
if($config{'general.ActiveDriver'}){
	print "Running ActiveDriver. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ACTIVEDRIVER/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/ACTIVEDRIVER/LATEST");

	# Generate config file
	generateConfig("ActiveDriver");

	# Run ActiveDriver
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_ActiveDriver -e $logsDir/ActiveDriver.error.log -o $logsDir/ActiveDriver.run.log $config{'general.scriptsDir'}/ACTIVEDRIVER/run_ActiveDriver.pl --config $analysisDir/ActiveDriver_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse ActiveDriver output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_ActiveDriver_parseOutput -e $logsDir/ActiveDriver_parseOutput.error.log -o $logsDir/ActiveDriver_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ACTIVEDRIVER/parse_to_standard_output.pl --inDir $analysisDir --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# FATHMM
if($config{'general.FATHMM'}){
	print "Running FATHMM. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/FATHMM/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/FATHMM/LATEST");

	# Generate config file
	generateConfig("FATHMM");

	# Run FATHMM
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_FATHMM -e $logsDir/FATHMM.error.log -o $logsDir/FATHMM.run.log $config{'general.scriptsDir'}/FATHMM/run_FATHMM.pl --config $analysisDir/FATHMM_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse FATHMM output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	my $numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_FATHMM_parseOutput -e $logsDir/FATHMM_parseOutput.error.log -o $logsDir/FATHMM_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/FATHMM/parse_to_standard_output.pl --in $analysisDir/fathmm.result --samples $numSamples --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# transFIC
if($config{'general.transFIC'}){
	print "Running transFIC. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/TRANSFIC/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/TRANSFIC/LATEST");

	# Generate config file
	generateConfig("transFIC");

	# Run transFIC
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_transFIC -e $logsDir/transFIC.error.log -o $logsDir/transFIC.run.log $config{'general.scriptsDir'}/TRANSFIC/run_transFIC.pl --config $analysisDir/transFIC_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse transFIC output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	my $numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_transFIC_parseOutput -e $logsDir/transFIC_parseOutput.error.log -o $logsDir/transFIC_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/TRANSFIC/parse_to_standard_output.pl --in $analysisDir/transFIC.result --samples $numSamples --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}


# HotNet2A
if($config{'general.HotNet2A'}){
	print "Running HotNet2A. Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/HOTNET2A/$runID";
	system("mkdir -p $analysisDir") unless (-e $analysisDir);
	system("ln -sfn $analysisDir $config{'general.analysisDir'}/HOTNET2A/LATEST");

	# Generate config file
	generateConfig("HotNet2A");

	# Run HotNet2
	if($mutSigID eq ""){
		$mutSigID = 9999;
	}
	$command = "$qsub -l mem_free=200G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_HotNet2A -e $logsDir/HotNet2A.error.log -o $logsDir/HotNet2A.run.log -hold_jid $mutSigID $config{'general.scriptsDir'}/HOTNET2A/run_HotNet2.pl --config $analysisDir/HotNet2A_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	# Parse HotNet2 output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_HotNet2A_parseOutput -e $logsDir/HotNet2A_parseOutput.error.log -o $logsDir/HotNet2A_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/HOTNET2A/parse_to_standard_output.pl --in $analysisDir/HotNet2.result --outDir $resultsDir";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}

# oncoIMPACT-discovery
if($config{'general.oncoIMPACT-discovery'}){
	print "Running oncoIMPACT (discovery mode). Please wait...";

	# Initialise folder
	$analysisDir = "$config{'general.analysisDir'}/ONCOIMPACT-DISCOVERY/$runID";
	unless (-e $analysisDir){
		system("mkdir -p $analysisDir");
		system("ln -sfn $analysisDir $config{'general.analysisDir'}/ONCOIMPACT-DISCOVERY/LATEST");
	}

	# Generate config file
	generateConfig("oncoIMPACT-discovery");

	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -hold_jid $oncoIMPACT_ID -pe 1 $config{'cluster.numThreads'} -N $config{'general.disease'}_oncoIMPACT-discovery -e $logsDir/oncoIMPACT-discovery.error.log -o $logsDir/oncoIMPACT-discovery.run.log $config{'general.scriptsDir'}/ONCOIMPACT/oncoIMPACT.pl $analysisDir/oncoIMPACT-discovery_$runID.cfg";
	$command = $command . " 1" if ($flag_debug);
	submit($command);

	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT-discovery_parseOutput -e $logsDir/oncoIMPACT-discovery_parseOutput.error.log -o $logsDir/oncoIMPACT-discovery_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/ONCOIMPACT/parse_to_standard_output.pl --in $analysisDir --out $resultsDir/oncoIMPACT-discovery.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);

	print "Job submitted.\n";
}




## Generate status report
$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_getStatus -e $logsDir/getStatus.error.log -o $logsDir/getStatus.run.log -hold_jid " . join (",", @queue) . " $config{'general.scriptsDir'}/get_status.pl --config $configFile";

close(TRACE);

### Sub-routines ###
sub submit {
	$command = "@_";
	print TRACE "[Command] $command\n";
	unless($flag_simulate){
		my $return = `$command`;
		chomp($return);
		push( @queue, $return );
		print TRACE "\tJob $return submitted.\n";
	}
}    # end sub submit

sub generateConfig {
	my $software = "@_";
	open(OUT, "> $analysisDir/${software}_$runID.cfg");
	print TRACE "[Config] Generating file: $analysisDir/${software}_$runID.cfg\n" if ($flag_debug);

	given($software){
		when( 'oncoIMPACT'){
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptDir=$config{'general.scriptsDir'}/ONCOIMPACT\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "cnv=$config{'oncoIMPACT.cnv'}\n";
			print OUT "exp=$config{'oncoIMPACT.exp'}\n";
			print OUT "snp=$config{'oncoIMPACT.snp'}\n";
			print OUT "dataType=$config{'oncoIMPACT.dataType'}\n";
			print OUT "databaseExport=$config{'oncoIMPACT.databaseExport'}\n";
			continue;
		}
		when( 'DriverNet' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "expData=$config{'DriverNet.expData'}\n";
			print OUT "mutMatrix=$config{'DriverNet.mutMatrix'}\n";
			print OUT "numProc=$config{'cluster.numThreads'}\n";
			print OUT "scriptsDir=$config{'DriverNet.tcgaDir'}\n";
			continue;
		}
		when( 'MutSigCV' ){
			my $temp = $analysisDir;
			print OUT "outDir=$temp\n";
			print OUT "matlab=$config{'MutSigCV.matlab'}\n";
			print OUT "maf=$config{'MutSigCV.maf'}\n";
			print OUT "coverage=$config{'MutSigCV.coverage'}\n";
			print OUT "covariate=$config{'MutSigCV.covariate'}\n";
			print OUT "dict=$config{'MutSigCV.dict'}\n";
			print OUT "chr=$config{'MutSigCV.chr'}\n";
			print OUT "prefix=$config{'general.disease'}\n";
			continue;
		}
		when( 'OncodriveFM' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'OncodriveFM.annotation'}\n";
			print OUT "mappingFile=$config{'OncodriveFM.mappingFile'}\n";
			print OUT "numThreads=1\n";
			continue;
		}
		when( 'DawnRank' ){
			print OUT "adj=$config{'DawnRank.adj'}\n";
			print OUT "exp=$config{'DawnRank.exp'}\n";
			print OUT "mut=$config{'DawnRank.mut'}\n";
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/DAWNRANK\n";
			continue;
		}
		when( 'OncodriveCLUST' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'OncodriveCLUST.annotation'}\n";
			print OUT "transcript=$config{'OncodriveCLUST.transcript'}\n";
			continue;
		}
		when( 'NetBox' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "gene=$config{'NetBox.gene'}\n";
			print OUT "mutationFrequency=$config{'NetBox.mutationFrequency'}\n";
			print OUT "maxMutation=$config{'NetBox.maxMutation'}\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/NETBOX\n";
			continue;
		}
		when( 'CHASM' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "maf=$config{'CHASM.maf'}\n";
			print OUT "hg18=$config{'CHASM.hg18'}\n";
			print OUT "classifier=$config{'CHASM.classifier'}\n";
			continue;
		}
		when( 'oncoIMPACT-v1'){
			print OUT "outDir=$analysisDir\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "cnv=$config{'oncoIMPACT-v1.cnv'}\n";
			print OUT "exp=$config{'oncoIMPACT-v1.exp'}\n";
			print OUT "snp=$config{'oncoIMPACT-v1.snp'}\n";
			print OUT "network=$config{'oncoIMPACT-v1.network'}\n";
			print OUT "benchmarkGeneList=$config{'oncoIMPACT-v1.benchmarkGeneList'}\n";
			print OUT "dataType=$config{'oncoIMPACT-v1.dataType'}\n";
			continue;
		}
		when( 'HotNet2' ){
			print OUT "mem=$config{'cluster.mem'}\n";
			print OUT "runtime=$config{'cluster.runtime'}\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "resultsDir=$analysisDir\n";
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/HOTNET2\n";
			print OUT "installationDir=$config{'HotNet2.installationDir'}\n";
			print OUT "sigThreshold=$config{'HotNet2.sigThreshold'}\n";
			print OUT "freqFile=$config{'HotNet2.freqFile'}\n";
			print OUT "irefMaf=$config{'HotNet2.irefMaf'}\n";
			print OUT "irefIndex=$config{'HotNet2.irefIndex'}\n";
			print OUT "irefNetworks=$config{'HotNet2.irefNetworks'}\n";
			print OUT "deltaPerm=$config{'HotNet2.deltaPerm'}\n";
			print OUT "sigPerm=$config{'HotNet2.sigPerm'}\n";
			continue;
		}
		when( 'S2N' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/S2N\n";
			print OUT "completeSamples=$config{'general.completeSamples'}\n";
			print OUT "exp=$config{'S2N.exp'}\n";
			print OUT "cnv=$config{'S2N.cnv'}\n";
			continue;
		}
		when( 'OncodriveCIS' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/ONCODRIVECIS\n";
			print OUT "disease=$config{'general.disease'}\n";
			print OUT "completeSamples=$config{'general.completeSamples'}\n";
			print OUT "exp=$config{'OncodriveCIS.exp'}\n";
			print OUT "cnv=$config{'OncodriveCIS.cnv'}\n";
			print OUT "normals=$config{'OncodriveCIS.normals'}\n";
			continue;
		}
		when( 'ActiveDriver' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'ActiveDriver.annotation'}\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/ACTIVEDRIVER\n";
			continue;
		}
		when( 'FATHMM' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'FATHMM.annotation'}\n";
			print OUT "ref=$config{'FATHMM.ref'}\n";
			continue;
		}
		when( 'transFIC' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'transFIC.annotation'}\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/TRANSFIC\n";
			continue;
		}
		when( 'HotNet2A' ){
			print OUT "mem=$config{'cluster.mem'}\n";
			print OUT "runtime=$config{'cluster.runtime'}\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "resultsDir=$analysisDir\n";
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptsDir=$config{'general.scriptsDir'}/HOTNET2\n";
			print OUT "installationDir=$config{'HotNet2A.installationDir'}\n";
			print OUT "sigThreshold=$config{'HotNet2A.sigThreshold'}\n";
			print OUT "freqFile=$config{'HotNet2A.freqFile'}\n";
			print OUT "mutSig=$config{'HotNet2A.mutSig'}\n";
			print OUT "irefMaf=$config{'HotNet2A.irefMaf'}\n";
			print OUT "irefIndex=$config{'HotNet2A.irefIndex'}\n";
			print OUT "irefNetworks=$config{'HotNet2A.irefNetworks'}\n";
			print OUT "deltaPerm=$config{'HotNet2A.deltaPerm'}\n";
			print OUT "sigPerm=$config{'HotNet2A.sigPerm'}\n";
			continue;
		}
		when( 'oncoIMPACT-discovery'){
			print OUT "outDir=$analysisDir\n";
			print OUT "scriptDir=$config{'general.scriptsDir'}/ONCOIMPACT\n";
			print OUT "numThreads=1\n";
			print OUT "cnv=$config{'oncoIMPACT-discovery.cnv'}\n";
			print OUT "exp=$config{'oncoIMPACT-discovery.exp'}\n";
			print OUT "snp=$config{'oncoIMPACT-discovery.snp'}\n";
			print OUT "dataType=$config{'oncoIMPACT-discovery.dataType'}\n";
			print OUT "database=$config{'oncoIMPACT-discovery.database'}\n";
			continue;
		}
		default{
			close(OUT);
		}
	}
}
