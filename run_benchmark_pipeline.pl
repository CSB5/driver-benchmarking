#!/usr/bin/perl

#use warnings;
use Config::Simple;
use Getopt::Long;
use POSIX 'strftime';
use 5.010;

my $version = "2.0.0";
my $date = strftime '%Y%m%d', localtime;
my $runID = "${version}_${date}";

my ( $configFile, $flag_debug, $flag_help, %config, @queue, $command, $lastID );
my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH";

my $help_message = "
This script runs various softwares for benchmarking.

Usage:
	run_benchmark_pipeline.pl [OPTIONS]

Options:
	--config = path to config file *
	--debug: prints trace to STDERR
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
system("mkdir -p $config{'general.analysisDir'}") unless (-d $config{'general.analysisDir'});
system("mkdir -p $config{'general.resultsDir'}") unless (-d $config{'general.resultsDir'});
print "done.\n";

# Initializing variables
print "Initializing variables. Please wait...";
system("mkdir -p $config{'general.analysisDir'}/LOGS") unless (-d "$config{'general.analysisDir'}/LOGS");
open( TRACE, ">$config{'general.analysisDir'}/LOGS/trace.log" );
@queue = ();
my $runtime = "$config{'cluster.runtime'}:0:0";
print "done.\n";



### Softwares/Packages ###
# oncoIMPACT
if($config{'general.oncoIMPACT'}){
	print "Running oncoIMPACT. Please wait...";
	
	# Run oncoIMPACT
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'oncoIMPACT.numThreads'} -N $config{'general.disease'}_oncoIMPACT -e $config{'general.analysisDir'}/LOGS/oncoIMPACT.error.log -o $config{'general.analysisDir'}/LOGS/oncoIMPACT.run.log $config{'oncoIMPACT.scriptsDir'}/oncoIMPACT.pl $config{'oncoIMPACT.configFile'} $config{'oncoIMPACT.subSample'}";
	$command = $command . " 1" if ($flag_debug);
	submit($command);
	
	# Parse oncoIMPACT output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT_parseOutput -e $config{'general.analysisDir'}/LOGS/oncoIMPACT_parseOutput.error.log -o $config{'general.analysisDir'}/LOGS/oncoIMPACT_parseOutput.run.log -hold_jid $lastID $config{'oncoIMPACT.scriptsDir'}/parse_to_standard_output.pl --in $config{'general.analysisDir'}/oncoIMPACT_analysis/GENE_LIST/ALTERATION.dat --out $config{'general.resultsDir'}/oncoIMPACT.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# DriverNet
if($config{'general.DriverNet'}){
	print "Running DriverNet. Please wait...";
	
	# Run DriverNet
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'DriverNet.numThreads'} -N $config{'general.disease'}_DriverNet -e $config{'general.analysisDir'}/LOGS/DriverNet.error.log -o $config{'general.analysisDir'}/LOGS/DriverNet.run.log $config{'DriverNet.scriptsDir'}/run_driver_net.pl --config $config{'DriverNet.configFile'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse DriverNet output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_DriverNet_parseOutput -e $config{'general.analysisDir'}/LOGS/DriverNet_parseOutput.error.log -o $config{'general.analysisDir'}/LOGS/DriverNet_parseOutput.run.log -hold_jid $lastID $config{'DriverNet.scriptsDir'}/parse_to_standard_output.pl --in $config{'general.analysisDir'}/DRIVERNET/res_driver_net.dat --out $config{'general.resultsDir'}/DriverNet.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}


# MutSigCV
if($config{'general.MutSigCV'}){
	print "Running MutSigCV. Please wait...";
	
	# Filter TCGA MAF file
	if($config{'MutSigCV.flagFilter'}){
		$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV_filterMAF -e $config{'general.analysisDir'}/LOGS/MutSigCV_filterMAF.error.log -o $config{'general.analysisDir'}/LOGS/MutSigCV_filterMAF.run.log $config{'DriverNet.scriptsDir'}/filter_maf.pl --samples $config{'general.completeSamples'} --maf $config{'MutSigCV.maf'} --outDir $config{'MutSigCV.dataDir'}";
		$command = $command . " --debug" if ($flag_debug);
		submit($command);
		
		$lastID = $queue[-1];
		chomp($lastID);
	} else{
		$lastID = 99999;
	}
	
	# Run MutSigCV
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV -e $config{'general.analysisDir'}/LOGS/MutSigCV.error.log -o $config{'general.analysisDir'}/LOGS/MutSigCV.run.log -hold_jid $lastID $config{'MutSigCV.scriptsDir'}/run_MutSigCV.pl --config $config{'MutSigCV.configFile'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	# Parse MutSigCV output to standard format
	$lastID = $queue[-1];
	chomp($lastID);
	$command = "$qsub -l mem_free=1G,h_rt=0:10:0 -pe OpenMP 1 -N $config{'general.disease'}_MutSigCV_parseOutput -e $config{'general.analysisDir'}/LOGS/MutSigCV_parseOutput.error.log -o $config{'general.analysisDir'}/LOGS/MutSigCV_parseOutput.run.log -hold_jid $lastID $config{'MutSigCV.scriptsDir'}/parse_to_standard_output.pl --in $config{'general.analysisDir'}/MUTSIGCV/$config{'general.disease'}.sig_genes.txt --out $config{'general.resultsDir'}/MutSigCV.result";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}



### Sub-routines ###
sub submit {
	$command = "@_";
	print TRACE "[Command] $command\n";
	my $return = `$command`;
	chomp($return);	
	push( @queue, $return );
	print TRACE "\tJob $return submitted.\n";
}    # end sub submit


sub generateConfig {
	my $software = "@_";
	open(OUT, "> $config{'general.analysisDir'}/$runID/${software}_$runID.cfg");
	
	given($software){
		when( 'oncoIMPACT'){
			print OUT "outDir=$config{'general.analysisDir'}/$runID\n";
			print OUT "scriptDir=$config{'oncoIMPACT.scriptsDir'}\n";
			print OUT "numThreads=$config{'cluster.numThreads'}\n";
			print OUT "cnv=$config{'oncoIMPACT.cnv'}\n";
			print OUT "exp=$config{'oncoIMPACT.exp'}\n";
			print OUT "snp=$config{'oncoIMPACT.snp'}\n";
			continue;
		}
		when( 'DriverNet' ){
			print OUT "dataDir=$config{'DriverNet.dataDir'}\n";
			print OUT "outDir=$config{'general.analysisDir'}/$runID/$software\n";
			print OUT "expData=$config{'DriverNet.expData'}\n";
			print OUT "numProc=$config{'cluster.numThreads'}\n";
			print OUT "scriptDir=$config{'DriverNet.tcgaDir'}\n";
			continue;
		}
		when( 'MutSigCV' ){
			print OUT "outDir=$config{'general.analysisDir'}/$runID/$software\n";
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