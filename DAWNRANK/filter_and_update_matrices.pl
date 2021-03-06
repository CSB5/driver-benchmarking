#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($expr_matrix_file, $mut_matrix_file, $network_file, $outDir, $flag_debug, $flag_help);

my $help_message = "
This script filters and updates the expression and mutation matrices required by DawnRank.

Usage:
	filter_and_update_matrices.pl [OPTIONS]

Options:
	--exp = path to expression matrix  *
	--mut = path to mutation matrix  *
	--adj = path to adjacency network matrix *
	--out = path to output directory *
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
	"exp=s"      	=> \$expr_matrix_file,
	"mut=s"      	=> \$mut_matrix_file,
	"adj=s"      	=> \$network_file,
	"out=s"			=> \$outDir,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "EXPRESSION: $expr_matrix_file\n";
	print STDERR "MUTATION: $mut_matrix_file\n";
	print STDERR "NETWORK: $network_file\n";
}

#Get all the gene in the network
my %gene_to_ID = ();
my @ID_to_gene = ();
my $last_ID = 0;

print "Reading network file. Please wait...";
open(FILE, $network_file);
<FILE>;#To skip the header
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    for(my $i = 0; $i < @line; $i++){
	$gene = $line[$i];
	if(! exists $gene_to_ID{$gene}){
	    push(@ID_to_gene, $gene);
	    $gene_to_ID{$gene} = $last_ID;
	    $last_ID++;
	}
    }
}
close(FILE);
print "done.\n";


#For expression data
#Contruct a patient x mutation gene matrix:
# * gene order similar as the network
# * all genes in the network present in the matrix

print "Reading expression matrix file. Pleae wait...";
open(FILE, "$expr_matrix_file");

my $first = 1;
my @sample_order = ();
my %sample_expression = ();
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    
    #Get the header with the sample names
    if($first){
	@sample_order = @line[1..@line-1];
	#print STDERR " **** @sample_order\n";
	foreach $s (@sample_order){
	    my %map = ();
	    $sample_expression{$s} = \%map;
	}
	$first = 0;
    }
    else{
	$gene = $line[0];
	next if(! exists $gene_to_ID{$gene});
	for(my $i = 0; $i < @sample_order; $i++){
	    $s = $sample_order[$i];
	    #if($line[$i + 1] == 1){
	    $sample_expression{$s}->{$gene} = abs($line[$i + 1]);
	    #}
	}
    }
}
close(FILE);
print "done.\n";


#Write the file
print "Generating updated expression matrix file. Please wait...";
open(OUT, "> $outDir/expression.dat");
print OUT "".join("\t", @sample_order)."\n";
my $str = "";
for(my $i = 0; $i < @ID_to_gene; $i++){
    $gene = $ID_to_gene[$i];
    $str = $gene;
    for(my $j = 0; $j < @sample_order; $j++){
	$s = $sample_order[$j];
	$v = 0;
	if(exists $sample_expression{$s}->{$gene}){
	    $v = $sample_expression{$s}->{$gene};
	}
	$str .= "\t".$v;
    }
    print OUT $str."\n";
}
close(OUT);
print "done.\n";


#For expression data
#Contruct a patient x mutation gene matrix:
# * gene order similar as the network
# * all genes in the network present in the matrix
# * same sample order as the one in expression matrix

#Contruct a patient x mutation gene matrix:
# * gene order similar as the network
# * all genes in the network present in the matrix

print "Reading mutation matrix file. Please wait...";
open(FILE, "$mut_matrix_file");

$first = 1;
my %sample_mutation = ();
my @gene_order = ();
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    
    #Get the header with the sample names
    if($first){
	@gene_order = @line;
	$first = 0;
    }
    else{
	$s = $line[0];
	for(my $i = 0; $i < @gene_order; $i++){
	    $gene = $gene_order[$i];
	    if($line[$i + 1] == 1){
		$sample_mutation{$s}->{$gene} = 1;
	    }
	}
    }
}
close(FILE);
print "done.\n";


#Write the file
print "Generating updated mutation matrix file. Please wait...";
open(OUT, "> $outDir/mutation.dat");
print OUT "".join("\t", @sample_order)."\n";
for(my $i = 0; $i < @ID_to_gene; $i++){
    $gene = $ID_to_gene[$i];
    $str = $gene;
    for(my $j = 0; $j < @sample_order; $j++){
	$s = $sample_order[$j];
	$v = 0;
	if(exists $sample_mutation{$s}->{$gene}){
	    $v =1;
	}
	$str .= "\t".$v;
    }
    print OUT $str."\n";
}
close(OUT);
print "done.\n";
