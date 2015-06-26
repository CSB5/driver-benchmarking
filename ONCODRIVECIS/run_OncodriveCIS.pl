#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs OncodriveCIS.

Usage:
	run_OncodriveCIS.pl [OPTIONS]

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

$command = "mkdir $config{'default.outDir'}";
print STDERR "$command\n" if($flag_debug);
system($command) unless (-d $config{'default.outDir'});

print "Preparing files needed for analysis. Please wait...";
#Read the completeSamples file to get the sample list
open(FILE, $config{'default.completeSamples'});
my $sample;
my %sample_list = ();
my @sample_order = ();
while(<FILE>){
    $sample = $_;
    chomp($sample);
    $sample_list{$sample} = 1;
}
close(FILE);


#read the gistic CNV file and write the cnv input FILE
print STDERR " *** Read the GIST2 CNV file $config{'default.cnv'}\n" if($flag_debug);
open(FILE, $config{'default.cnv'});
#Get the GISTIC sample order
$header = <FILE>;chop $header;
@tab = split(/\t/, $header);
for(my $i = 3; $i < @tab; $i++){
    $full_sample_name = $tab[$i];
    #Need to convert the sample name;
    $sample_name = convert_sample_name($full_sample_name);
    if(exists $sample_list{$sample_name}){
		$sample_list{$sample_name} = $i;
    }
}

#Read the matrix and write the cnv input FILE
#-2 and 2 CNV region => 1 CNA region
#-1 or 1 CNV +1 => 0 uncertain
#0 or not present diploid region => 0 diploid
print STDERR " *** Output the CNV file $config{'default.outDir'}/cnv.gistic.gene.tsv\n" if($flag_debug);
open(OUT, ">$config{'default.outDir'}/cnv.gistic.gene.tsv");
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene = $line[0];
    foreach $s (keys %sample_list){
		$cnv_val = $line[$sample_list{$s}];
		if($cnv_val != 0){
		    print OUT $gene."\t".$s."\t".$cnv_val."\n";
		}
    }
}
close(FILE);
close(OUT);

#Read the normal sample file and output the sample to use file
print STDERR " *** Out put the sample used file\n";
open(FILE, $config{'default.normals'});
while(<FILE>){
    $full_sample_name = $_;chop $full_sample_name;
    #Convert the name
    $sample_name = convert_sample_name($full_sample_name);
    $sample_list{$sample_name} = -1;
}
open(OUT, ">$config{'default.outDir'}/samples_to_process.tsv");
foreach $s (keys %sample_list){
    $val = 1;
    $val = 0 if($sample_list{$s} == -1);
    print OUT $s."\t".$val."\n";
}

#Use of the deSEQ matrix that contains all the normals
print STDERR " *** Read Expression File $config{'default.exp'}\n" if($flag_debug);
open(FILE, $config{'default.exp'});
$header = <FILE>;chop $header;
@tab = split(/\t/, $header);
for(my $i = 0; $i < @tab; $i++){
    $full_sample_name = $tab[$i];
    #Need to convert the sample name;
    $sample_name = convert_sample_name($full_sample_name);
    #print STDERR " *** the sample $i -> $sample_name\n";<STDIN>;
    if(exists $sample_list{$sample_name}){
	$sample_list{$sample_name} = $i;
	push(@sample_order, $sample_name);
    }
}
print STDERR " *** Output the expression file\n";
open(OUT, ">$config{'default.outDir'}/expression.gene.tsv");
#Write the header
$header = "gene";
for(my $i = 0; $i < @sample_order; $i++){
    $sample = $sample_order[$i];
    if(exists $sample_list{$sample}){
	$header .= "\t".$sample;
    }
}
print OUT $header."\n";

my $log_2 = log(2);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene = $line[0];
    print OUT $gene;
    for(my $i = 0; $i < @sample_order; $i++){
	$sample = $sample_order[$i];
	#print STDERR " *** the sample $i -> $sample\n";<STDIN>;
	if(exists $sample_list{$sample}){
	    $val = log($line[$i+1]+0.001)/$log_2;
	    print OUT "\t".$val;
	}
    }
    print OUT "\n";
}
close(FILE);
close(OUT);

print "done.\n";


# Run OncodriveCIS
print "Running OncodriveCIS. Please wait...";
$command = "oncodrivecis.py -e $config{'default.outDir'}/expression.gene.tsv -c $config{'default.outDir'}/cnv.gistic.gene.tsv -s $config{'default.outDir'}/samples_to_process.tsv -o $config{'default.outDir'} -p $config{'default.disease'}";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";


# Compute pVal & FDR
print "Computing pVal & FDR. Please wait...";
$command = "$config{'default.scriptsDir'}/compute_pVal_fdr.R $config{'default.outDir'}";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";


#### Sub-routines
sub convert_sample_name{
    my ($f_s_n) = @_;
    @tmp = split(/\-/, $f_s_n);
    $s_n = join("-", @tmp[0..3]);
    return $s_n;
}
