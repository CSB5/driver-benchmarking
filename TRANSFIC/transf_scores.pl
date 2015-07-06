#!/usr/bin/perl

use strict 'vars';
use File::Basename;

#########################################################################################

my $partition = $ARGV[0];
my $mutsfile = $ARGV[1];
my $outfile = $ARGV[2];

die "Usage:\n\t./transf_scores.pl\t<partition (gosmf|cp|gosbp|doms)>\t<Mutations input file>\t[output file]\n\n\tMutations file format\n\t\tMutation id (optional)\tGene ID (Must be Ensembl id)\tSIFT score\tPolyphen2 score\tMutationAssessor score\n\n" if (scalar @ARGV < 2);

#########################################################################################

my $comm_check = "head -1 $mutsfile";
my $check = `$comm_check`;
my @fields = split /\s+/, $check;

if (scalar @fields < 4) {
  die "Incorrect mutations file format!\n\n";
}
else {
  if ((scalar @fields == 4 && $fields[0] !~ /ENS[G|T|P]\d+/) || (scalar @fields == 5 && $fields[1] !~ /ENS[G|T|P]\d+/)){
    die "Please convert your identifiers to Ensembl ids\n\n";
  }
}

open (MUTS, $mutsfile);
my @muts = <MUTS>;

if ($outfile) {
  open (STDOUT, ">$outfile");
}

my $dist_file = dirname($0) . "/data/$partition.genes";

my (%original, %transformed);

for (my $i = 0; $i < @muts; $i++){
  my @data = split /\s+/, $muts[$i];

  if (scalar @fields == 4){
    $original{gene} = $data[0];
    $original{sift} = $data[1];
    $original{pph2} = $data[2];
    $original{ma} = $data[3];
    $original{id}++;
  }
  elsif (scalar @fields == 5){
    $original{gene} = $data[1];
    $original{sift} = $data[2];
    $original{pph2} = $data[3];
    $original{ma} = $data[4];
    $original{id} = $data[0];
  }

  my $get_lines = "grep $original{gene} $dist_file | cut -f2,3,4,5,6,7";
  my $line = `$get_lines`;
  chomp $line;

  if ($line !~ /\d+/){
    print STDERR "Unable to transform scores of mutation $original{id}; $original{gene} not in precomputed database of scores distributions...\n";
    next;
  }

  my %genes2dist;

  my @dists = split /\t/, $line;
  $genes2dist{siftmean} = $dists[0];
  $genes2dist{siftstd} = $dists[1];
  $genes2dist{pph2mean} = $dists[2];
  $genes2dist{pph2std} = $dists[3];
  $genes2dist{mamean} = $dists[4];
  $genes2dist{mastd} = $dists[5];

  $transformed{id} = $original{id};
  $transformed{gene} = $original{gene};

  if ($original{sift} == 0){
    $original{sift} = 0.001;
  }
  elsif ($original{sift} == 1){
    $original{sift} = 0.999;
  }

  $original{sift} = log(1 - $original{sift})-log($original{sift});
  $transformed{sift} = sprintf ("%.3f", ($original{sift} - $genes2dist{siftmean})/$genes2dist{siftstd});

  if ($original{pph2} == 0){
    $original{pph2} = 0.001;
  }
  elsif ($original{pph2} == 1){
    $original{pph2} = 0.999;
  }

  $original{pph2} = log($original{pph2})-log(1 - $original{pph2});
  $transformed{pph2} = sprintf ("%.3f", ($original{pph2} - $genes2dist{pph2mean})/$genes2dist{pph2std});

  $transformed{ma} = sprintf ("%.3f", ($original{ma} - $genes2dist{mamean})/$genes2dist{mastd});

  my %class;

  if ($transformed{sift} < -1){
    $class{sift} = 'low_impact';
  }
  elsif ($transformed{sift} >= -1 && $transformed{sift} < 2){
    $class{sift} = 'medium_impact';
  }
  else {
    $class{sift} = 'high_impact';
  }

  if ($transformed{pph2} < -1){
    $class{pph2} = 'low_impact';
  }
  elsif ($transformed{pph2} >= -1 && $transformed{pph2} < 1.5){
    $class{pph2} = 'medium_impact';
  }
  else {
    $class{pph2} = 'high_impact';
  }

  if ($transformed{ma} < -1){
    $class{ma} = 'low_impact';
  }
  elsif ($transformed{ma} >= -1 && $transformed{ma} < 2){
    $class{ma} = 'medium_impact';
  }
  else {
    $class{ma} = 'high_impact';
  }

  print "$transformed{id}\t$transformed{gene}\t$transformed{sift}\t$class{sift}\t$transformed{pph2}\t$class{pph2}\t$transformed{ma}\t$class{ma}\n";
}
