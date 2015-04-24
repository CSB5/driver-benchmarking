#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;
use JSON;
use HTTP::Request::Common;
use LWP::UserAgent;
use File::Fetch;
use Data::Dumper;


my ($configFile, %config, $flag_debug, $flag_help);

my $help_message = "
This script runs CHASM via the CRAVAT web API.

Usage:
	run_CHASM.pl [OPTIONS]

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

my $url = "http://www.cravat.us/rest/service";
my $email = "chiakhb\@gis.a-star.edu.sg";
my $waitTime = 1800; # time to wait between checking status of jobs (in seconds) [default = 30 minutes]

my ($jobID, $inputFile, $fileURL, $status);
$inputFile = "$config{'default.outDir'}/chasm_input.dat";

# Generate input file
print "Generating CHASM compatible input file. Please wait...";
generateInput();
print "done.\n";

# Submit job
print "Submitting job. Please wait...";
$jobID = submitJob();
print "done. JobID:$jobID\n";

# Check status
print "Checking status of job. Current job status:\n";
while($status eq "Running"){
	sleep $waitTime;
	checkJob();
}

# Get results
print "Downloading results. Please wait...";
getResults();
print "done.\n";


#### Sub-routines ####
sub generateInput {
	my (@temp, $index);
	open(IN, $config{'default.maf'});
	open(OUT, ">$inputFile");
	<IN>;	# skip header
	$index = 1;
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		if( $temp[4] =~ /^chr/ ){
			print OUT "$index\t$temp[4]\t$temp[5]\t$temp[7]\t$temp[10]\t$temp[12]\t$temp[15]\n";
		} else{
			print OUT "$index\tchr$temp[4]\t$temp[5]\t$temp[7]\t$temp[10]\t$temp[12]\t$temp[15]\n";
		}
		$index++;
	}
	close(IN);
	close(OUT);
} #end generateInput

sub submitJob {
	my $ua = LWP::UserAgent->new;
	my $req = $ua->request(	POST "$url/submit",
	              			Content_Type => 'form-data',
	              			Content => [
								analyses => 'CHASM',
								analysistype => 'driver',
								chasmclassifier => $config{'default.classifier'},
								email => $email,
								functionalannotation => 'off',
								hg18 => 'off',
								inputfile => ["$inputFile"],
								mupitinput => 'off',
								tsvreport => 'on'
							]
						);
	
	if($flag_debug){
		print STDERR "## Job submission ##\n";
		print STDERR Dumper($req);
		print STDERR "\n";
	}
	
	# Check the outcome of the response
	if ($req->is_success) {
		my $data = decode_json($req->content);
		$status = "Running";
		return $data->{'jobid'};
	}
	else {
		die "\n NOT SUCCESSFUL\n";
	}	
} #end submitJob

sub checkJob {
	my $ua = LWP::UserAgent->new;
	my $req = $ua->request(GET "$url/status?jobid=$jobID");
	my $data = decode_json($req->content);
	
	if($flag_debug){
		print STDERR "## Job status ##\n";
		print STDERR Dumper($req);
		print STDERR "\n";
	}
	
	if($data->{'status'} eq "Success" ){
		$fileURL = $data->{'resultfileurl'};
		$status = "Success";
		print "\tCompleted.\n";
	} else{
		$status = "Running";
		print "\tRunning.\n";
	}	
} #end checkJob

sub getResults {
	my $ff = File::Fetch->new(uri => $fileURL );
	my $where = $ff->fetch( to => $config{'default.outDir'} ) or die $ff->error;
	system("unzip -j -d $config{'default.outDir'} $config{'default.outDir'}/$jobID.zip");
	system("rm -rf $config{'default.outDir'}/$jobID.zip");	
} #end getResults