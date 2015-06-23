#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ($configFile, %config, $flag_debug, $flag_help, $command);

my $help_message = "
This script prepares the necessary inputs and runs S2N.

Usage:
	run_S2N.pl [OPTIONS]

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


#Get the list used for the analysis
print "Preparing files needed for analysis. Please wait...";
my $gene_list_file = "$config{'default.outDir'}/gene_list_order.txt";
my $gene_list_expression_file = "$config{'default.outDir'}/gene_list_expression.txt";

$command = "tail -n+2 $config{'default.exp'} | cut -f 1 > $gene_list_expression_file";
print STDERR "$command\n" if($flag_debug);
system($command);

$command = "cut -f 1 $gene_list_expression_file $config{'default.cnv'} | grep -v \"Gene\" | sort | uniq -d > $gene_list_file";
print STDERR "$command\n" if($flag_debug);
system($command);

my @gene_order = ();
my %gene_list = ();
open(FILE, $gene_list_file);
while(<FILE>){
    $gene = $_;chop $gene;
    $gene_list{$gene} = 1;
}
close(FILE);

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
    #print STDERR " *** $sample_name\n";<STDIN>;
    if(exists $sample_list{$sample_name}){
		$sample_list{$sample_name} = $i;
		push(@sample_order, $sample_name);
		#print STDERR " *** $sample_name\t$i\n";<STDIN>;
    }
}

#Read the matrix and write the cnv input FILE
#-2 and 2 CNV region => 1 CNA region
#-1 or 1 CNV +1 => 0 uncertain
#0 or not present diploid region => 0 diploid
print STDERR " *** Output the CNV file $config{'default.outDir'}/cnv.gistic.gene.tsv\n" if($flag_debug);
open(OUT, ">$config{'default.outDir'}/cnv.gistic.gene.tsv");
print OUT "".(join("\t", @sample_order))."\n";
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene = $line[0];
    next if(! exists $gene_list{$gene});
    push(@gene_order, $gene);
    print OUT $gene;
    for(my $i = 0; $i < @sample_order; $i++){
		$s = $sample_order[$i];
		$cnv_val = $line[$sample_list{$s}];
		#print STDERR $gene."\t".$cnv_val."\t".$sample_list{$s}."\t".$s."\n";<STDIN>;
		$cnv_val = ($cnv_val == -2 || $cnv_val == 2) ? 1 : 0;
		print OUT "\t".$cnv_val;
    }
    print OUT "\n";
}
close(FILE);
close(OUT);

#sort the expression file according to gene list
print STDERR " *** Read Initial Expression File $config{'default.exp'}\n" if($flag_debug);
open(FILE, $config{'default.exp'});

#Get the expression header
my %gene_expression = ();
$header = <FILE>;chop $header;
my @expr_sample_order = split(/\t/, $header);

while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene = $line[0];
    for(my $i = 0; $i < @expr_sample_order; $i++){
    	$sample = $expr_sample_order[$i];
		next if(! exists $gene_list{$gene});
		#print STDERR " **** |$gene|\n";<STDIN>;
		$exp_val = $line[$i+1];
		if(! exists $gene_expression{$gene}){
		    my %map = ();
		    $gene_expression{$gene} = \%map;
		}
		$gene_expression{$gene}->{$sample} = $exp_val;
    }
}
close(FILE);


print STDERR " *** Output Expression File $config{'default.outDir'}/expression_reodered.dat\n" if($flag_debug);
open(OUT, ">$config{'default.outDir'}/expression_reodered.dat");
print OUT "".(join("\t", @sample_order))."\n";
for (my $i = 0; $i < @gene_order; $i++){
    $gene = $gene_order[$i];
    print OUT $gene;
    for(my $j = 0; $j < @sample_order; $j++){
		$sample = $sample_order[$j];
		print OUT "\t".$gene_expression{$gene}->{$sample};
    }
}

close(OUT);

print "done.\n";

# Run S2N in R
print "Running S2N. Please wait...";
$command = "$config{'default.scriptsDir'}/S2N.R $config{'default.outDir'}/expression_reodered.dat $config{'default.outDir'}/cnv.gistic.gene.tsv $config{'default.outDir'}";
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
