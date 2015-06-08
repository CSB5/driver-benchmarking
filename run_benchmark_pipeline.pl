#!/usr/bin/perl

#use warnings;
#no warnings 'experimental';
use Config::Simple;
use Getopt::Long;
use POSIX 'strftime';
use 5.010;

my $version = "v3.5.0";
my $date = strftime '%Y%m%d', localtime;
my $runID = "${date}_${version}";

my ( $configFile, $flag_debug, $flag_help, $flag_update, $flag_simulate, %config, @queue, $command, $lastID, $software, $outDir);
my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH,NETBOX_HOME";

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
	
	# CHASM
	print "CHASM: ";
	$command = "$config{'general.scriptsDir'}/CHASM/parse_to_standard_output.pl --inDir $config{'general.analysisDir'}/CHASM/LATEST --out $resultsDir/CHASM.result";
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
	
	# LJB
	print "LJB: ";
	my $numSamples = `wc -l $config{'general.completeSamples'} | cut -f 1 -d \" \"`;
	chomp($numSamples);
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
	
	# OncoIMPACT-v1
	print "OncoIMPACT-v1: ";
	$command = "$config{'general.scriptsDir'}/ONCOIMPACT-V1/parse_to_standard_output.pl --in $config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST/ --out $resultsDir/oncoIMPACT-v1.result";
	if(-s "$config{'general.analysisDir'}/ONCOIMPACT-V1/LATEST/driver_list.txt"){
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
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=23:0:0,h=n070 -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV -e $logsDir/MutSigCV.error.log -o $logsDir/MutSigCV.run.log $config{'general.scriptsDir'}/MUTSIGCV/run_MutSigCV.pl --config $analysisDir/MutSigCV_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse MutSigCV output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
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
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_OncodriveFM -e $logsDir/OncodriveFM.error.log -o $logsDir/OncodriveFM.run.log $config{'general.scriptsDir'}/ONCODRIVEFM/run_OncodriveFM.pl --config $analysisDir/OncodriveFM_$runID.cfg";
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
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_OncodriveCLUST -e $logsDir/OncodriveCLUST.error.log -o $logsDir/OncodriveCLUST.run.log $config{'general.scriptsDir'}/ONCODRIVECLUST/run_OncodriveCLUST.pl --config $analysisDir/OncodriveCLUST_$runID.cfg";
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
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_DawnRank -e $logsDir/DawnRank.error.log -o $logsDir/DawnRank.run.log $config{'general.scriptsDir'}/DAWNRANK/run_DawnRank.pl --config $analysisDir/DawnRank_$runID.cfg";
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
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_CHASM_parseOutput -e $logsDir/CHASM_parseOutput.error.log -o $logsDir/CHASM_parseOutput.run.log -hold_jid $lastID $config{'general.scriptsDir'}/CHASM/parse_to_standard_output.pl --inDir $analysisDir --out $resultsDir/CHASM.result";
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
			#print OUT "databaseExport=$config{'oncoIMPACT.databaseExport'}\n";
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
			print OUT "prefix=$config{'MutSigCV.prefix'}\n";
			continue;
		}		
		when( 'OncodriveFM' ){
			print OUT "outDir=$analysisDir\n";
			print OUT "annotation=$config{'OncodriveFM.annotation'}\n";
			print OUT "mappingFile=$config{'OncodriveFM.mappingFile'}\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
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
		default{
			close(OUT);
		}
	}
}