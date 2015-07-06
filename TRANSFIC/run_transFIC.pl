#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs transFIC.

Usage:
	run_transFIC.pl [OPTIONS]

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
my %score_col = (
  "SIFT_score", 12,
  "Polyphen2_HDIV_score", 14,
  "Polyphen2_HVAR_score", 16,
  "LRT_score", 18,
  "MutationTaster_score", 20,
  "MutationAssessor_score", 22,
  "FATHMM_score", 24,
  "RadialSVM_score", 26,
  "LR_score", 28,
  "VEST3_score", 30,
);

my %score_used = (
  "SIFT_score", "SIFT",
  "Polyphen2_HVAR_score", "PPH2",
  "MutationAssessor_score", "MA"
);
my @score_order = keys %score_used;
my ($file_in, $file_intermediate, $file_final);
$file_in = "$config{'default.outDir'}/transFIC.input";
$file_intermediate = "$config{'default.outDir'}/transFIC.intermediate";
$file_final = "$config{'default.outDir'}/transFIC.result";

open(OUT, ">$res_file");
open(FILE, $config{'default.annotation'});
while(<FILE>){
  chop $_;
  @line = split(/\t/, $_);
  next if($line[0] eq "Chr");
  $all_define = 1;
  my @res_score = ();
  for($i = 0; $i < @score_order; $i++){
	  $s = $score_order[$i];
	  $s_val = $line[$score_col{$s}];
	  push(@res_score, $s_val);
	  if($s_val eq "."){
      $all_define = 0;
	  }
  }

  if($all_define){
  	$sample_ID = $line[@line-1];

  	#Get the ensble gene ID
  	#$gene_ENS = $line[43];
  	@gene_ENS_list = split(",", $line[43]);
  	$gene_ENS = $gene_ENS_list[0];
  	#print STDERR " list size ".(@gene_ENS_list+0)." -> @gene_ENS_list -> $gene_ENS\n";<STDIN>;

  	$gene = $line[6];

  	print OUT $sample_ID."_".$gene."\t".$gene_ENS."\t".(join("\t", @res_score))."\n";
  }
}
close(FILE);
close(OUT);

print "done.\n";


# Run transFIC
print "Running transFIC. Please wait...";
$command = "$config{'default.scriptsDir'}/transf_scores.pl gosmf $file_in $file_intermediate";
print STDERR "$command\n" if($flag_debug);
system($command);

$command = "cut -f 1,7 $file_intermediate | sed 's/_/\\t/' > $file_final";
print STDERR "$command\n" if($flag_debug);
system($command);
print "done.\n";
