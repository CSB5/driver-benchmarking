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
	1.1

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
#There is now 2 reference files
$ref_file_ensembl = $config{'default.refEnsembl'};
$ref_file_uniprot = $config{'default.refUniprot'};
#
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
    # "stopgain", "stop",
    # "stoploss", "stop",
    );

#Notice that not all gene have a protein ID
#Gene without protein ID will be dropped for the analysis
#All trancrips with a valid protein will be given as imput to fathmm
#
#Read the ensembl file
my %trans_to_prot = ();
open(FILE, $ref_file_ensembl);
<FILE>; # skip header
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $trans = $line[0];
    
    #For the ensembl file
    $prot = $line[9];
    $gene = $line[10];
    
    chop $gene if(index($gene, ",") != -1);
    next if(index($gene, ",") != -1 || $gene eq "n/a");

    if($prot ne "" && index($prot, "NP_") == -1){
	my @tab = ($prot, $gene);
	$trans_to_prot{$trans} = \@tab;
    }
}
close(FILE);
#
#Read the uniprot file
open(FILE, $ref_file_uniprot);
<FILE>; # skip header
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $trans = $line[0];
    
    #For the ucs file
    $prot = $line[8];
    $gene = $line[9];
    
    chop $gene if(index($gene, ",") != -1);
    next if(index($gene, ",") != -1 || $gene eq "n/a");

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
my (@all_effect, @all_annot);
my ( $counter_pass, $counter_fail, $counter_snv, $counter_flag );
$counter_fail = 0;
$counter_snv = 0;
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);

    next if($line[0] eq "Chr");

    @splicing_effect = split(/\;/, $line[5]);
    @splicing_effect = ($line[5]) if(@splicing_effect == 0);
    
    @sample_full_name = split(/\-/, $line[@line-1]);
    $sample_name = join("-", @sample_full_name[0..3]);
    
    #For usc annotation
    @all_effect_ucsc = split(/\;/, $line[8]);
    #@all_effect_ucsc = ($line[8]) if(@all_effect == 0);
    $info_line = $line[9];#this line contain the gene coordinate data
    #
    @all_annot_ucsc = split(/\,/, $info_line);
    @all_annot_ucsc = ($info_line) if(@all_annot_ucsc == 0);
    #
    
    #For ensembl annotation
    @all_effect_ensembl = split(/\;/, $line[45]);
    @all_effect_ensembl = ($line[45]) if(@all_effect == 0);
    $info_line = $line[46];#this line contain the gene coordinate data
    #
    @all_annot_ensembl = split(/\,/, $info_line);
    @all_annot_ensembl = ($info_line) if(@all_annot_ensembl == 0);
    
    

    
    for($i = 0; $i < 2; $i++){
	#Select the annotation type
	@all_effect = @all_effect_ensembl;@all_annot = @all_annot_ensembl;
	if($i == 1){@all_effect = @all_effect_ucsc;@all_annot = @all_annot_ucsc;}
	#
	#fathmm do not take into account splicing event
	if($all_effect[0] eq "." || $all_effect[0] eq "unknown" || $splicing_effect[0] eq "splicing"){
	    next;
	}
	
	$counter_flag = 1;
	$counter_pass = 0;
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
		
		$counter_pass++;
		$counter_snv++ if( $counter_flag );
		$counter_flag = 0;
		
		$coord = $annot_info[4];
		#print STDERR "the coord: ".$coord."\n";
		if($coord =~ /p\.(.*)/){
		    if( index( $1, ";" ) != -1){
			$prot_coord = (split( /;/, $1 ))[0];
		    } else{
			$prot_coord = $1;
		    }
		    #print STDERR $_."\n".$prot_coord."\n";<STDIN>;
		}
		$gene_name = $trans_to_prot{$trans_name}->[1];
		$prot_name = $trans_to_prot{$trans_name}->[0];
		
		$res = $prot_name."\t".$prot_coord."\t".$gene_name."\t".$sample_name."\t".$prot_name.":".$prot_coord."\t".$line[0].":".$line[1]."\n";
		
		print OUT $res;
		
	    }
	}
    }
    $counter_fail++ if( $counter_pass == 0 && $counter_flag == 0);
}
close(FILE);
close(OUT);
print "done.\n";
print "## Statistics of transcript to protein mapping ##\n";
print "Number of SNV without any mapping: $counter_fail (" . ( $counter_fail / $counter_snv ) . ")\n";
$counter_pass = $counter_snv - $counter_fail;
print "Number of SNV with mapping: $counter_pass (" . ( $counter_pass / $counter_snv ) . ")\n";

# Running FATHMM
print "Running FATHMM. Please wait...";
$command = "/mnt/software/unstowable/fathmm/cgi-bin/fathmm.py -w Cancer $input_file $output_file";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";

#exit 0;

# generate gene prediction score
#print "Generating gene prediction score. Please wait...";
#$command = "awk '{print \$2\":\"\$3\"\\t\"\$4\"\\t\"\$5}' $output_file | join -1 1 -2 5 --nocheck-order - $input_file | sed 's/ /\\t/g' | awk '{if(\$2 == \"CANCER\") print \$6\"\\t\"\$2\"\\t\"\$3\"\\t\"\$7}' | sort -k1,1 -k4,4 > $result_file";
#print STDERR "$command\n" if(1 || $flag_debug);
#system($command);

#This was too messy and need to re-coded in perl

#Read the input file
open(FILE, $input_file);
my %mutation_ID_info;
my %snv_list_without_pred;
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    my @info = ($line[2], $line[3], $line[5]);#(gene, sample, snv)
    $ID = $line[4];

    #print STDERR $ID."\t".$line[5]."\n";<STDIN>;

    $mutation_ID_info{$ID} = \@info;
    $snv_list_without_pred{$line[5]} = 1;
}
close(FILE);

my %snv_list_with_pred;
#Read the fathmm result file
open(FILE, $output_file);
open(OUT, ">$result_file");
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    my $ID = $line[2].":".$line[3];
    my $prediction = $line[4];

    if($ID eq "ENSP00000341429:I106T"){
	#print STDERR " *** ".$ID."\t".$prediction."\n";<STDIN>;
    }

    if($prediction eq "CANCER" || $prediction eq "PASSENGER/OTHER"){
	if($prediction eq "CANCER"){
	    $score = $line[5];
	    $out_line = $mutation_ID_info{$ID}->[0]."\t".$prediction."\t".$score."\t".$mutation_ID_info{$ID}->[1]."\n";
	    print OUT $out_line;
	    if($ID eq "ENSP00000341429:I106T"){
		#print STDERR " *** ".$out_line."\n";<STDIN>;
	    }
	}
	$snv_ID = $mutation_ID_info{$ID}->[2];
	$snv_list_with_pred{$snv_ID} = 1;
	delete $snv_list_without_pred{$snv_ID};
    }
}
close(FILE);
close(OUT);

my $nb_snv_with_pred = keys(%snv_list_with_pred);
print " SNV with prediction: $nb_snv_with_pred ".($nb_snv_with_pred/$counter_snv)."\n";

my $nb_to_print =0;
foreach $g (keys %snv_list_without_pred){
    print STDERR $g."\n";
    $nb_to_print++;
    last if($nb_to_print == 10);
}


print "done.\n";

#ucsc ref
#/home/chiakhb/SCRIPTS/MUTATION_BENCHMARK/ONCODRIVECLUST/DATABASE/hg18_knownGene_withGeneName_withUniProt.txt
#ensembl ref
#/home/chiakhb/SCRIPTS/MUTATION_BENCHMARK/ONCODRIVECLUST/DATABASE/hg18_ensembl_withGeneName_withProteinID.txt
