#!/usr/bin/perl

#use warnings;
#no warnings 'experimental';
use Config::Simple;
use Getopt::Long;
use POSIX 'strftime';
use 5.010;

my $version = "2.0.2";
my $date = strftime '%Y%m%d', localtime;
my $runID = "${version}_${date}";

my ( $configFile, $flag_debug, $flag_help, $flag_simulate, %config, @queue, $command, $lastID, $software, $outDir );
my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH";

my $help_message = "
This script runs various softwares for benchmarking.

Usage:
	run_benchmark_pipeline.pl [OPTIONS]

Options:
	--config = path to config file *
	--debug: prints trace to STDERR
	--sim: simulation mode - only prints commands to be executed to TRACE log
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

# Preparing output folders
print "Preparing output folders. Please wait...";
my $analysisDir = "$config{'general.analysisDir'}/$runID";
my $resultsDir = "$analysisDir/CONSOLIDATED_RESULTS";
my $logsDir = "$analysisDir/LOGS";
system("mkdir -p $analysisDir") unless (-d $analysisDir);
system("mkdir -p $resultsDir") unless (-d $resultsDir);
system("mkdir -p $logsDir") unless (-d $logsDir);
print "done.\n";

# Initializing variables
print "Initializing variables. Please wait...";
open( TRACE, "> $analysisDir/LOGS/trace.log" );
@queue = ();
my $runtime = "$config{'cluster.runtime'}:0:0";
print "done.\n";



### Softwares/Packages ###
# oncoIMPACT
if($config{'general.oncoIMPACT'}){	
	print "Running oncoIMPACT. Please wait...";
	
	# Generate config file
	generateConfig("oncoIMPACT");
	
	
	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'oncoIMPACT.numThreads'} -N $config{'general.disease'}_oncoIMPACT -e $logsDir/oncoIMPACT.error.log -o $logsDir/oncoIMPACT.run.log $config{'oncoIMPACT.scriptsDir'}/oncoIMPACT.pl $analysisDir/oncoIMPACT_$runID.cfg $config{'oncoIMPACT.subSample'}";
	$command = $command . " 1" if ($flag_debug);
	submit($command);
	
	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT_parseOutput -e $logsDir/oncoIMPACT_parseOutput.error.log -o $logsDir/oncoIMPACT_parseOutput.run.log -hold_jid $lastID $config{'oncoIMPACT.scriptsDir'}/parse_to_standard_output.pl --in $analysisDir/oncoIMPACT_analysis/GENE_LIST/ALTERATION.dat --out $resultsDir/oncoIMPACT.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# DriverNet
if($config{'general.DriverNet'}){
	print "Running DriverNet. Please wait...";
	
	# Initialise folder
	$outDir = "$analysisDir/DRIVERNET";
	system("mkdir -p $outDir") unless (-d $outDir);
	
	# Generate config file
	generateConfig("DriverNet");
	
	# Run DriverNet
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'DriverNet.numThreads'} -N $config{'general.disease'}_DriverNet -e $logsDir/DriverNet.error.log -o $logsDir/DriverNet.run.log $config{'DriverNet.scriptsDir'}/run_driver_net.pl --config $analysisDir/DriverNet_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse DriverNet output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_DriverNet_parseOutput -e $logsDir/DriverNet_parseOutput.error.log -o $logsDir/DriverNet_parseOutput.run.log -hold_jid $lastID $config{'DriverNet.scriptsDir'}/parse_to_standard_output.pl --in $outDir/res_driver_net.dat --out $resultsDir/DriverNet.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# MutSigCV
if($config{'general.MutSigCV'}){
	print "Running MutSigCV. Please wait...";
	
	# Initialise folder	
	$outDir = "$analysisDir/MUTSIGCV";
	system("mkdir -p $outDir") unless (-d $outDir);
	
	# Generate config file
	generateConfig("MutSigCV");
	
	# Filter TCGA MAF file
	if($config{'MutSigCV.flagFilter'}){
		$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV_filterMAF -e $logsDir/MutSigCV_filterMAF.error.log -o $logsDir/MutSigCV_filterMAF.run.log $config{'MutSigCV.scriptsDir'}/filter_maf.pl --samples $config{'general.completeSamples'} --maf $config{'MutSigCV.maf'} --outDir $config{'MutSigCV.dataDir'}";
		$command = $command . " --debug" if ($flag_debug);
		submit($command);
		
		$lastID = $queue[-1];
		chomp($lastID);
	} else{
		$lastID = 99999;
	}
	
	# Run MutSigCV
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV -e $logsDir/MutSigCV.error.log -o $logsDir/MutSigCV.run.log -hold_jid $lastID $config{'MutSigCV.scriptsDir'}/run_MutSigCV.pl --config $analysisDir/MutSigCV_$runID.cfg'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse MutSigCV output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV_parseOutput -e $logsDir/MutSigCV_parseOutput.error.log -o $logsDir/MutSigCV_parseOutput.run.log -hold_jid $lastID $config{'MutSigCV.scriptsDir'}/parse_to_standard_output.pl --in $outDir/$config{'general.disease'}.sig_genes.txt --out $resultsDir/MutSigCV.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


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
			print OUT "scriptDir=$config{'oncoIMPACT.scriptsDir'}\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "cnv=$config{'oncoIMPACT.cnv'}\n";
			print OUT "exp=$config{'oncoIMPACT.exp'}\n";
			print OUT "snp=$config{'oncoIMPACT.snp'}\n";
			continue;
		}
		when( 'DriverNet' ){
			print OUT "dataDir=$config{'DriverNet.dataDir'}\n";
			print OUT "outDir=$analysisDir/$software\n";
			print OUT "expData=$config{'DriverNet.expData'}\n";
			print OUT "numProc=$config{'cluster.numThreads'}\n";
			print OUT "scriptDir=$config{'DriverNet.tcgaDir'}\n";
			continue;
		}
		when( 'MutSigCV' ){
			print OUT "outDir=$analysisDir/$software\n";
			print OUT "matlab=$config{'MutSigCV.matlab'}\n";
			print OUT "maf=$config{'MutSigCV.maf'}\n";
			print OUT "coverage=$config{'MutSigCV.coverage'}\n";
			print OUT "covariate=$config{'MutSigCV.covariate'}\n";
			print OUT "dict=$config{'MutSigCV.dict'}\n";
			print OUT "chr=$config{'MutSigCV.chr'}\n";
			print OUT "prefix=$config{'MutSigCV.prefix'}\n";
			continue;
		}
		default{
			close(OUT);
		}
	}
}