#!/usr/bin/perl
use warnings;

my ($result_file, $out_file) = @ARGV;

open(FILE, "sort -k2,2nr $result_file |");
open(OUT, ">$out_file");
print OUT "Gene_name\tSample\tRank\tpValue\tInfo\n";	# print header
my $rank = 1;
while(<FILE>){
    @line = split(/\t/, $_);
    $gene =  $line[0];
    $mut_freq = $line[1];
    print OUT $gene."\t"."ALL"."\t".$rank."\t".$mut_freq."\t"."-"."\n";
    $rank++;
}
close(FILE);
