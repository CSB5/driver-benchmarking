#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs FATHMM.

Usage:
	run_FATHMM.pl [OPTIONS]

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

# Generate FATHMM input files
my ($annot_file, $ref_file, $out_dir, $input_file, $output_file, $result_file);
$annot_file = $config{'default.annotation'};
$ref_file = $config{'default.ref'};
$out_dir = $config{'default.outDir'};
$input_file = "$out_dir/fathmm.input";
$output_file = "$out_dir/fathmm.output";
$result_file = "$out_dir/fathmm.result";

unless (-d $config{'default.outDir'}){
  $command = "mkdir $config{'default.outDir'}";
  print STDERR "$command\n" if($flag_debug);
  system($command);
}

print "Preparing files needed for analysis. Please wait...";
my %effect_name = (
    "nonsynonymous SNV", "non-synonymous",
    "stopgain", "stop",
    "stoploss", "stop",
    );

#Notice that not all gene have a protein ID
#Gene without protein ID will be dropped for the analysis
#All trancrips with a valid protein will be given as imput to fathmm
my %trans_to_prot = ();
open(FILE, $ref_file);
<FILE>; # skip header
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $trans = $line[0];
    $prot = $line[8];
    $gene = $line[9];

    if($prot ne "" && index($prot, "NP_") == -1){
	my @tab = ($prot, $gene);
	$trans_to_prot{$trans} = \@tab;
    }
}
close(FILE);


#Read the annovar file
open(FILE, $annot_file);
open(OUT, ">$input_file");
my $prot_coord = 0;
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);

    next if($line[0] eq "Chr");

    @splicing_effect = split(/\;/, $line[5]);
    @splicing_effect = ($line[5]) if(@splicing_effect == 0);

    @all_effect = split(/\;/, $line[8]);
    @all_effect = ($line[8]) if(@all_effect == 0);

    $info_line = $line[9];#this line contain the gene coordinate data

    #fathmm do not take into account splicing event
    if($all_effect[0] eq "." || $all_effect[0] eq "unknown" || $splicing_effect[0] eq "splicing"){
	next;
    }


    @sample_full_name = split(/\-/, $line[@line-1]);
    $sample_name = join("-", @sample_full_name[0..3]);

    @all_annot = split(/\,/, $info_line);
    @all_annot = ($info_line) if(@all_annot == 0);
    #

    for(my $i = 0; $i < @all_annot; $i++){
	$annot = $all_annot[$i];
	@annot_info = split(/\:/, $annot);
	$trans_name = $annot_info[1];
	if(!defined $trans_name){
	    print STDERR " *** Warning no transcript name available $all_effect[0]\n".$_;#<STDIN>;
	}
	#print STDERR $trans_name."\n";
	if(exists $trans_to_prot{$trans_name}){
	    #The effect representatin is quite weird
	    $eff = $all_effect[0];#We sre actually only looking at the first effect repported that should be the most important one
	    #$eff = $all_effect[$i] if(@all_effect != 1);

	    $effect = $effect_name{$eff};
	    if(! defined $effect){
		next;
	    }

	    $coord = $annot_info[4];
	    #print STDERR "the coord: ".$coord."\n";
	    if($coord =~ /p\.(.*)/){
		$prot_coord = $1;
		#print STDERR $_."\n".$prot_coord."\n";<STDIN>;
	    }
	    $gene_name = $trans_to_prot{$trans_name}->[1];
	    $prot_name = $trans_to_prot{$trans_name}->[0];

	    $res = $prot_name."\t".$prot_coord."\t".$gene_name."\t".$sample_name."\t".$prot_name.":".$prot_coord."\n";

	    print OUT $res;

	}
    }
}
close(FILE);
close(OUT);
print "done.\n";

# Running FATHMM
print "Running FATHMM. Please wait...";
$command = "/mnt/software/unstowable/fathmm/cgi-bin/fathmm.py -w Cancer $input_file $output_file";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";

# generate gene prediction score
print "Generating gene prediction score. Please wait...";
$command = "awk '{print \$2\":\"\$3\"\\t\"\$4\"\\t\"\$5}' $output_file | join -1 1 -2 5 --nocheck-order - $input_file | sed 's/ /\\t/g' | awk '{if(\$2 == \"CANCER\") print \$6\"\\t\"\$2\"\\t\"\$3\"\\t\"\$7}' > $result_file";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";
