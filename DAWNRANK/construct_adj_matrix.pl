#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($edge_list_file, $out_file, $flag_debug, $flag_help);

my $help_message = "
This script constructs the adjacency matrix required by DawnRank.

Usage:
	construct_adj_matrix.pl [OPTIONS]

Options:
	--adj = path to adjacency network matrix file *
	--out = path to output file *
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
	"adj=s"      	=> \$edge_list_file,
	"out=s"         => \$out_file,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "EDGE LIST: $edge_list_file\n";
	print STDERR "OUTPUT: $out_file\n";	
}


#my $nb_pair = 0;
#my $max_pair = 20;

my %gene_to_ID = ();
my @ID_to_gene = ();
my @adj_matrix = ();
my $last_ID = 0;

print "Constructing adjacency matrix. Please wait...";

open(FILE, $edge_list_file);

while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    
    #Update the IDs and the adjacency matrix
    for($i = 0; $i < @line; $i++){
	$gene = $line[$i];

	if(! exists $gene_to_ID{$gene}){
	    #The ID
	    push(@ID_to_gene, $gene);
	    $gene_to_ID{$gene} = $last_ID;

	    #The adjacency matrix
	    my @tab = ();
	    push(@adj_matrix, \@tab);
	    
	    $last_ID++;
	    
	}
    }
    
    #Add the edge to the adjancy matrix
    $index1 = $gene_to_ID{$line[0]};
    $index2 = $gene_to_ID{$line[1]};

    $adj_matrix[$index1]->[$index2] = 1;
    $adj_matrix[$index2]->[$index1] = 1;
    
    #$nb_pair++;
    #last if($nb_pair == $max_pair);

}

#output 
open(OUT, ">$out_file");
#The header with the gene names
print OUT "".(join("\t", @ID_to_gene))."\n";

my $str = "";
for(my $i = 0; $i < $last_ID; $i++){
    $str = $ID_to_gene[$i];
    for(my $j = 0; $j < $last_ID; $j++){
	$v = $adj_matrix[$i]->[$j];
	if(! defined $v){
	    $v = 0;
	}
	$str .= "\t".$v;
    }
    print OUT $str."\n";
}

print "done\n";


