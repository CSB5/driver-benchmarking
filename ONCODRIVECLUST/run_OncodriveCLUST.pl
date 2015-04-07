#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script runs OncodriveCLUST.

Usage:
	run_OncodriveCLUST.pl [OPTIONS]

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


# Generate OncodriveCLUST input files
my ($annot_file, $transcript_file, $out_dir);
$annot_file = $config{'default.annotation'};
$transcript_file = $config{'default.transcript'};
$out_dir = $config{'default.outDir'}; 

my %effect_name = (
    "nonsynonymous SNV", "non-synonymous",
    "stopgain", "stop",
    "stoploss", "stop",
    "splicing", "splice",
    "synonymous SNV", "synonymous"
    );

my %all_trascript = ();
open(FILE, $transcript_file);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene_name = $line[0];
    $trans_name = $line[1];
    $all_trascript{$trans_name} = $gene_name;
}
close(FILE);

$header = "symbol\tgene\ttranscript\tSample\tct\tposition";
open(OUT_C, ">$out_dir/nonsyn.txt");
open(OUT_S, ">$out_dir/syn.txt");
print OUT_C $header."\n";
#print OUT_S $header."\n";

my $nb_snv = 0;
my $nb_un_annotated_snv = 0;
my $nb_snv_without_effect = 0;
my $nb_snv_without_transcript = 0;
my $info_line;
my $splicing_offset;#the annotation of splicing variant miss the gene name -_-'
my @all_effect;

open(FILE, $annot_file);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    
    next if($line[0] eq "Chr");
    
    $nb_snv++;
    
    @splicing_effect = split(/\;/, $line[5]);
    @splicing_effect = ($line[5]) if(@all_effect == 0);

    @all_effect = split(/\;/, $line[8]);
    @all_effect = ($line[8]) if(@all_effect == 0);
    $splicing_offset = 0;

    $info_line = $line[9];#this line contain the gene coordinate data

    if($all_effect[0] eq "." && $splicing_effect[0] eq "splicing"){
	$info_line = $line[7];
	@all_effect = @splicing_effect;
	$splicing_offset = 1;
    }
    	
    
    @sample_full_name = split(/\-/, $line[@line-1]);
    $sample_name = join("-", @sample_full_name[0..3]);

    @all_annot = split(/\,/, $info_line);
    @all_annot = ($info_line) if(@all_annot == 0);
    #
    
    #

    #This could be due to splicing event, need to take into account this thing
    if($all_annot[0] eq "." || $all_annot[0] eq "UNKNOWN"){
	$nb_un_annotated_snv++;
	next;
    }

    $flag_t = 0;
    $flag_e = 0;
    for(my $i = 0; $i < @all_annot; $i++){
	$annot = $all_annot[$i];
	@annot_info = split(/\:/, $annot);
	$trans_name = $annot_info[1 - $splicing_offset];
	#print STDERR $trans_name."\n";
	if(exists $all_trascript{$trans_name}){
	    $flag_t = 1;
	    #The effect representatin is quite weird
	    $eff = $all_effect[0];#We sre actually only looking at the first effect repported that should be the most important one
	    #$eff = $all_effect[$i] if(@all_effect != 1);
	    $effect = $effect_name{$eff};
	    #
	    if(! defined $effect){
		next;
	    }
	    $flag_e = 1;
	    $coord = $annot_info[3 - $splicing_offset];
	    #print STDERR "the coord: ".$coord."\n";
	    if($coord =~ /(.*)(A|T|G|C)(\d+)(A|T|G|C)/){
		$cds_coord = $3;
		#print STDERR $_."\n".$cds_coord."\n";<STDIN>;
	    }
	    $gene_name = $all_trascript{$trans_name};
	    
	    $res = $gene_name."\t".$gene_name."\t".$trans_name."\t".$sample_name."\t".$effect."\t".$cds_coord."\n";

	    if($effect eq "synonymous"){
		print OUT_S $res;
	    }
	    else{
		print OUT_C $res;
	    }
	    last;
	}
    }
    if($flag_t == 0){
	$nb_snv_without_transcript++;
	#print STDERR " *** SNV without transcript:\n$_\n *** INFO LINE|$info_line|\n";<STDIN>;
    }
    else{
	if($flag_e == 0){
	    $nb_snv_without_effect++;
	}
    }

}

close(OUT_C);
close(OUT_S);

print STDERR " *** nb SNV analysised: $nb_snv:"
    ."\n\t".($nb_un_annotated_snv/$nb_snv)." non annoted SNVs"
    ."\n\t".($nb_snv_without_transcript/$nb_snv)." SNVs without transcript"
    ."\n\t".($nb_snv_without_effect/$nb_snv)." SNVs without effect"."\n";

close(FILE);


# Run OncodriveCLUST
$command = "export PATH=/mnt/software/unstowable/anaconda/bin/:\$PATH; source activate oncodriveclust; oncodriveclust -m 3 -o $out_dir/ $out_dir/nonsyn.txt $out_dir/syn.txt $transcript_file";
print STDERR "$command\n" if $flag_debug;
system($command);

