#!/usr/bin/perl

use warnings;
no warnings 'experimental';
use Config::Simple;
use Getopt::Long;
use POSIX 'strftime';
use 5.010;

my $version = "2.0.0";
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
$software = "oncoIMPACT"; 
if($config{'general.$software'}){
	print "Running $software. Please wait...";

	# Generate config file
	generateConfig($software);
	
	
	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'$software.numThreads'} -N $config{'general.disease'}_${software} -e $logsDir/$software.error.log -o $logsDir/$software.run.log $config{'$software.scriptsDir'}/oncoIMPACT.pl $analysisDir/${software}_$runID.cfg $config{'$software.subSample'}";
	$command = $command . " 1" if ($flag_debug);
	submit($command);
	
	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_${software}_parseOutput -e $logsDir/${software}_parseOutput.error.log -o $logsDir/${software}_parseOutput.run.log -hold_jid $lastID $config{'$software.scriptsDir'}/parse_to_standard_output.pl --in $analysisDir/oncoIMPACT_analysis/GENE_LIST/ALTERATION.dat --out $resultsDir/$software.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# DriverNet
$software = "DriverNet"; 
if($config{'general.$software'}){
	print "Running $software. Please wait...";
	
	# Initialise folder
	$outDir = "$analysisDir/DRIVERNET";
	system("mkdir -p $outDir") unless (-d $outDir);
	
	# Generate config file
	generateConfig($software);
	
	# Run DriverNet
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'$software.numThreads'} -N $config{'general.disease'}_$software -e $logsDir/$software.error.log -o $logsDir/$software.run.log $config{'$software.scriptsDir'}/run_driver_net.pl --config $analysisDir/${software}_$runID.cfg";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse DriverNet output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_${software}_parseOutput -e $logsDir/${software}_parseOutput.error.log -o $logsDir/${software}_parseOutput.run.log -hold_jid $lastID $config{'$software.scriptsDir'}/parse_to_standard_output.pl --in $outDir/res_driver_net.dat --out $resultsDir/$software.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# MutSigCV
$software = "MutSigCV"; 
if($config{'general.$software'}){
	print "Running $software. Please wait...";
	
	# Initialise folder	
	$outDir = "$analysisDir/MUTSIGCV";
	system("mkdir -p $outDir") unless (-d $outDir);
	
	# Generate config file
	generateConfig($software);
	
	# Filter TCGA MAF file
	if($config{'$software.flagFilter'}){
		$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_${software}_filterMAF -e $logsDir/${software}_filterMAF.error.log -o $logsDir/${software}_filterMAF.run.log $config{'$software.scriptsDir'}/filter_maf.pl --samples $config{'general.completeSamples'} --maf $config{'$software.maf'} --outDir $config{'$software.dataDir'}";
		$command = $command . " --debug" if ($flag_debug);
		submit($command);
		
		$lastID = $queue[-1];
		chomp($lastID);
	} else{
		$lastID = 99999;
	}
	
	# Run MutSigCV
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_$software -e $logsDir/$software.error.log -o $logsDir/$software.run.log -hold_jid $lastID $config{'$software.scriptsDir'}/run_MutSigCV.pl --config $analysisDir/${software}_$runID.cfg'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse MutSigCV output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_${software}_parseOutput -e $logsDir/${software}_parseOutput.error.log -o $logsDir/${software}_parseOutput.run.log -hold_jid $lastID $config{'$software.scriptsDir'}/parse_to_standard_output.pl --in $outDir/$config{'general.disease'}.sig_genes.txt --out $resultsDir/$software.result";
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