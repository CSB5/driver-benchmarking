#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::Util qw(sum);

my ($file_in, $numSamples, $outDir, $file_out, $flag_debug, $flag_help);

my $help_message = "
This script parses ANNOVAR's annotation output to a standard output for LJB methods.

Usage:
	parse_to_standard_output.pl [OPTIONS]

Options:
	--in = path to ANNOVAR's annotation results file *
	--samples = number of samples *
	--outDir = path to output folder *
	--debug: prints trace to STDERR
	--help : prints this message

* indicates required parameters


Version:
	1.2

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}

GetOptions(
	"in=s"      	=> \$file_in,
	"samples=i"		=> \$numSamples,
	"outDir=s"      => \$outDir,
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


my ($counter, @temp, $gene, $score, $type, $threshold, %reportedGenes, $currentGene, @currentScore, $currentSample, @samples, $sample, $bestScore);

sift("7,13,58");
polyphen2("7,17,58");
mutationAssessor("7,23,24,58");
mutationTaster("7,21,22,58");


######## sub routines ########
# SIFT
sub sift {
	$threshold = 0.05;
	$file_out = "$outDir/SIFT.temp";
	open(OUT, ">$file_out");
	open(IN, "cut -f $_[0] $file_in | tail -n+2 | sort -k1,1 -k3,3 |");
	$currentGene = "";
	@currentScore = ();
	$currentSample = "";
	$bestScore = 1;
	@samples = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$sample = $temp[2];
		next if ($score eq "." || $score > $threshold);	# skip current entry as <score> doesn't meet requirement
		if( $gene ne $currentGene && $currentGene ne ""){	# new gene encountered
			# push previous entry to array
			push(@currentScore, (1 - $bestScore));
			push(@samples, $currentSample);

			# Print results and re-initialize variables
			print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t" . join(";", @samples) . "\n";
			$currentGene = "";
			@currentScore = ();
			$currentSample = "";
			$bestScore = $score;
			@samples = ();
		}
		$currentSample = $sample if ($currentSample eq "");
		$currentGene = $gene if ($currentGene eq "");
		if( $sample ne $currentSample){	# new sample encountered
			# push previous sample to array and re-initialize variables
			push(@currentScore, (1 - $bestScore));
			push(@samples, $currentSample);
			$currentSample = $sample;
			$bestScore = 1;
		}
		$bestScore = $score if ($score < $bestScore);
	}
	# Print last entry
	push(@currentScore, (1 - $bestScore));
	push(@samples, $currentSample);
	print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t" . join(";", @samples) . "\n";
	close(OUT);
	close(IN);

	$counter = 1;
	$file_out = "$outDir/SIFT.result";
	open(IN, "sort -k2,2gr $outDir/SIFT.temp |");
	open(OUT, ">$file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$sample = $temp[2];
		print OUT $gene . "\t" . $sample . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/SIFT.temp") unless $flag_debug;
}


# Polyphen2_HVAR_score
sub polyphen2 {
	# $threshold = 0.909; # stringent mode
	$threshold = 0.446; # relaxed mode (possibly damaging)
	$file_out = "$outDir/PolyPhen2.temp";
	open(OUT, ">$file_out");
	open(IN, "cut -f $_[0] $file_in | tail -n+2 | sort -k1,1 -k3,3 |");
	$currentGene = "";
	@currentScore = ();
	$currentSample = "";
	$bestScore = 0;
	@samples = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$sample = $temp[2];
		next if ($score eq "." || $score < $threshold);	# skip current entry as <score> doesn't meet requirement
		if( $gene ne $currentGene && $currentGene ne ""){	# new gene encountered
			# push previous entry to array
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);

			# Print results and re-initialize variables
			print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t" . join(";", @samples) . "\n";
			$currentGene = "";
			@currentScore = ();
			$currentSample = "";
			$bestScore = $score;
			@samples = ();
		}
		$currentSample = $sample if ($currentSample eq "");
		$currentGene = $gene if ($currentGene eq "");
		if( $sample ne $currentSample){	# new sample encountered
			# push previous sample to array and re-initialize variables
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);
			$currentSample = $sample;
			$bestScore = 0;
		}
		$bestScore = $score if ($score > $bestScore);
	}
	# Print last entry
	push(@currentScore, $bestScore);
	push(@samples, $currentSample);
	print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t" . join(";", @samples) . "\n";
	close(OUT);
	close(IN);

	$counter = 1;
	$file_out = "$outDir/PolyPhen2.result";
	open(IN, "sort -k2,2gr $outDir/PolyPhen2.temp |");
	open(OUT, ">$file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$sample = $temp[2];
		print OUT $gene . "\t" . $sample . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/PolyPhen2.temp") unless $flag_debug;
}


# MutationAssessor
sub mutationAssessor {
	$file_out = "$outDir/MutationAssessor.temp";
	open(OUT, ">$file_out");
	open(IN, "cut -f $_[0] $file_in | tail -n+2 | sort -k1,1 -k4,4 |");
	$currentGene = "";
	@currentScore = ();
	$currentSample = "";
	$bestScore = 0;
	@samples = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		$sample = $temp[3];
		next if ($score eq ".");	# skip current entry as <score> doesn't meet requirement
		next unless ($type eq "H"	# stringent mode (high only)
					 || $type eq "M" # relaxed mode (high and medium only)
		);	# skip current entry as <type> doesn't meet requirement
		if( $gene ne $currentGene && $currentGene ne ""){	# new gene encountered
			# push previous entry to array
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);

			# Print results and re-initialize variables
			print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t". $type .  "\t" . join(";", @samples) . "\n";
			$currentGene = "";
			@currentScore = ();
			$currentSample = "";
			$bestScore = $score;
			@samples = ();
		}
		$currentSample = $sample if ($currentSample eq "");
		$currentGene = $gene if ($currentGene eq "");
		if( $sample ne $currentSample){	# new sample encountered
			# push previous sample to array and re-initialize variables
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);
			$currentSample = $sample;
			$bestScore = 0;
		}
		$bestScore = $score if ($score > $bestScore);
	}
	# Print last entry
	push(@currentScore, $bestScore);
	push(@samples, $currentSample);
	print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t". $type .  "\t" . join(";", @samples) . "\n";
	close(OUT);
	close(IN);

	$counter = 1;
	$file_out = "$outDir/MutationAssessor.result";
	open(IN, "sort -k2,2gr $outDir/MutationAssessor.temp | ");
	open(OUT, ">$file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		$sample = $temp[3];
		next if ($score eq ".");
		print OUT $gene . "\t" . $sample . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$reportedGenes{$gene} = "";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/MutationAssessor.temp") unless $flag_debug;
}



# MutationTaster
sub mutationTaster {
	$file_out = "$outDir/MutationTaster.temp";
	open(OUT, ">$file_out");
	open(IN, "cut -f $_[0] $file_in | tail -n+2 | sort -k1,1 -k4,4 |");
	$currentGene = "";
	@currentScore = ();
	$currentSample = "";
	$bestScore = 0;
	@samples = ();
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		$sample = $temp[3];
		next if ($score eq ".");	# skip current entry as <score> doesn't meet requirement
		next unless ($type eq "A" || $type eq "D");	# skip current entry as <type> doesn't meet requirement
		if( $gene ne $currentGene && $currentGene ne ""){	# new gene encountered
			# push previous entry to array
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);

			# Print results and re-initialize variables
			print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t". $type .  "\t" . join(";", @samples) . "\n";
			$currentGene = "";
			@currentScore = ();
			$currentSample = "";
			$bestScore = $score;
			@samples = ();
		}
		$currentSample = $sample if ($currentSample eq "");
		$currentGene = $gene if ($currentGene eq "");
		if( $sample ne $currentSample){	# new sample encountered
			# push previous sample to array and re-initialize variables
			push(@currentScore, $bestScore);
			push(@samples, $currentSample);
			$currentSample = $sample;
			$bestScore = 0;
		}
		$bestScore = $score if ($score > $bestScore);
	}
	# Print last entry
	push(@currentScore, $bestScore);
	push(@samples, $currentSample);
	print OUT $currentGene . "\t" . sum(@currentScore)/$numSamples . "\t". $type .  "\t" . join(";", @samples) . "\n";
	close(OUT);
	close(IN);

	$counter = 1;
	$file_out = "$outDir/MutationTaster.result";
	open(IN, "sort -k2,2gr $outDir/MutationTaster.temp |");
	open(OUT, ">$file_out");
	print OUT "Gene_name\tSample\tRank\tScore\tInfo\n";	# print header
	while(<IN>){
		chomp(@temp = split(/\t/, $_));
		$gene = $temp[0];
		$score = $temp[1];
		$type = $temp[2];
		$sample = $temp[3];
		print OUT $gene . "\t" . $sample . "\t" . $counter . "\t" . $score . "\t" . "-" . "\n";
		$reportedGenes{$gene} = "";
		$counter++;
	}
	close(OUT);
	close(IN);
	system("rm -f $outDir/MutationTaster.temp") unless $flag_debug;
}
