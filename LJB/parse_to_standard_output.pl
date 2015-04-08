#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::Util qw(sum);

my ($file_in, $outDir, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses ANNOVAR's annotation output to a standard output for LJB methods.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to ANNOVAR's annotation results file *
	--outDir = path to output folder *
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
	"in=s"      	=> \$file_in,
	"outDir=s"         => \$outDir,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "INPUT: $file_in\n";
	print STDERR "OUTPUT DIR: $outDir\n";	
}


my ($counter, @temp, $gene, $score, $type, $threshold, %reportedGenes, $currentGene, @currentScore);

sift();
polyphen2();
mutationAssessor();
mutationTaster();

######## sub routines ########
sub mean {
    return sum(@_)/@_;
}

# SIFT
sub sift {	
	$threshold = 0.05;
	$file_out = "$outDir/SIFT.temp";
	open(OUT, "> $file_out");
	open(IN, "cut -f 7,13 $file_in | tail -n+2 | sort -k1,1 |");
	$currentGene = "";
	@currentScore = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		next if ($score eq "." || $score > $threshold);
		if( $gene ne $currentGene && $currentGene ne ""){
			print OUT $gene . "\t" . mean(@currentScore) . "\n";
			$currentGene = "";
			@currentScore = ();				
		} else{
			$currentGene = $gene if ($currentGene eq "");
			push(@currentScore, $score);
		}		
	}
	close(OUT);
	close(IN);
		
	$counter = 1;
	$file_out = "$outDir/SIFT.result";
	open(IN, "sort -k2,2n $outDir/SIFT.temp |");
	open(OUT, "> $file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/SIFT.temp");
}


# Polyphen2_HVAR_score
sub polyphen2 {
	$threshold = 0.909; # stringent mode
	# $threshold = 0.446; # relaxed mode
	$file_out = "$outDir/PolyPhen2.temp";
	open(OUT, "> $file_out");
	open(IN, "cut -f 7,17 $file_in | tail -n+2 | sort -k1,1 |");
	$currentGene = "";
	@currentScore = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		next if ($score eq "." || $score < $threshold);
		if( $gene ne $currentGene && $currentGene ne ""){
			print OUT $gene . "\t" . mean(@currentScore) . "\n";
			$currentGene = "";
			@currentScore = ();				
		} else{
			$currentGene = $gene if ($currentGene eq "");
			push(@currentScore, $score);
		}		
	}
	close(OUT);
	close(IN);
	
	$counter = 1;
	$file_out = "$outDir/PolyPhen2.result";
	open(IN, "sort -k2,2nr $outDir/PolyPhen2.temp |");
	open(OUT, "> $file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/PolyPhen2.temp");
}


# MutationAssessor
sub mutationAssessor {
	$file_out = "$outDir/MutationAssessor.temp";
	open(OUT, "> $file_out");
	open(IN, "cut -f 7,23,24 $file_in | tail -n+2 | sort -k1,1 |");
	$currentGene = "";
	@currentScore = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];		
		$type = $temp[2];
		next if ($score eq ".");
		next unless ($type eq "H"	# stringent mode
					 #|| $type eq "M" # relaxed mode
		);
		if( $gene ne $currentGene && $currentGene ne ""){
			print OUT $gene . "\t" . mean(@currentScore) . "\t". $type . "\n";
			$currentGene = "";
			@currentScore = ();				
		} else{
			$currentGene = $gene if ($currentGene eq "");
			push(@currentScore, $score);
		}		
	}
	close(OUT);
	close(IN);
	
	$counter = 1;
	$file_out = "$outDir/MutationAssessor.result";
	open(IN, "sort -k2,2nr $outDir/MutationAssessor.temp | ");
	open(OUT, "> $file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		next if ($score eq ".");
		print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$reportedGenes{$gene} = "";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/MutationAssessor.temp");
}


# MutationTaster
sub mutationTaster {
	$file_out = "$outDir/MutationTaster.temp";
	open(OUT, "> $file_out");
	open(IN, "cut -f 7,21,22 $file_in | tail -n+2 | sort -k1,1 |");
	$currentGene = "";
	@currentScore = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];		
		$type = $temp[2];
		next if ($score eq ".");
		next unless ($type eq "A" || $type eq "D");
		if( $gene ne $currentGene && $currentGene ne ""){
			print OUT $gene . "\t" . mean(@currentScore) . "\t". $type . "\n";
			$currentGene = "";
			@currentScore = ();				
		} else{
			$currentGene = $gene if ($currentGene eq "");
			push(@currentScore, $score);
		}		
	}
	close(OUT);
	close(IN);
	
	$counter = 1;
	$file_out = "$outDir/MutationTaster.result";
	open(IN, "sort -k2,2nr $outDir/MutationTaster.temp |");
	open(OUT, "> $file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		print OUT $gene . "\t" . "ALL" . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$reportedGenes{$gene} = "";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/MutationTaster.temp");
}