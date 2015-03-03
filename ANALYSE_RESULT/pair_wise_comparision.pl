#!/usr/bin/perl
use warnings;

my ($consolidated_result_dir, $data_gene_annotation_file, $out_dir, $gene_status_selection) = @ARGV;

`mkdir $out_dir`;

$gene_status_selection = "" if(! defined $gene_status_selection);

my %threshold_list = 
    ("MUT_FREQ", [0.10, 0.05, 0.02, 0.01, 0],
     "RANK", [5, 10, 20, 50, 100]);
     #"MUT_FREQ",  [0]);
     #"RANK", [5]);

my %method_ID = ();
my @ID_to_method = ();
my %sample_list = ();

#
#To update using that file:
#/mnt/pnsg10_projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/TEST_DATA/COAD/gene_mutation_frequency.txt
#

my (%pan_cancer, %cancer_census);
cancer_annotation();

open(FILE, $data_gene_annotation_file);
%data_gene_annotation = ();
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $gene =  $line[0];
    
    #Gene annotation
    $mut_freq = $line[1];

    $gene_status = "-";
    $gene_status = "CANCER" if(exists $pan_cancer{$gene} || exists $cancer_census{$gene});

    #Sample annotation
    @sample_info_tmp = split(/\;/, $line[2]);
    my @sample_info = ();
    foreach $s_m (@sample_info_tmp){
	@tmp = split(/\:/, $s_m);
	$sample = $tmp[0];
	$sample_list{$sample} = 1;
	push(@sample_info, $sample);
    }
    
    #Update the geneannotation
    my %gene_annot = ("MUT_FREQ", $mut_freq, "STATUS", $gene_status, "SAMPLE", \@sample_info);
    $data_gene_annotation{$gene} = \%gene_annot;
    
}
close(FILE);

opendir(DIR, $consolidated_result_dir);
@all_result_file = readdir(DIR);
close(DIR);

my %method_result = ();

#Construct the data base of results
foreach $res_file (@all_result_file){
    if($res_file =~ m/(.*)\.result/){
	#The method analysed
	$method = $1;
	$method_ID{$method} = @ID_to_method+0;
	push(@ID_to_method, $method);
	
	#To get the info for the method
	my %method_info = ();
	open(FILE, "$consolidated_result_dir/$res_file");
	<FILE>;#Skip the header
	while(<FILE>){
	    chop $_;
	    @line = split(/\t/, $_);
	    $gene = $line[0];
	    $rank = $line[2];
	    my %gene_info = ("RANK", $rank, "FREQ_CALL", 0, "USE", 0);
	    $method_info{$gene} = \%gene_info;
	}
	close(FILE);
	$method_result{$method} = \%method_info;
    }
}

#Init the pairwise matrise method distance
my @matrix = ();
for(my $i = 0; $i < @ID_to_method; $i++){
    my @tab = ();
    for(my $j = 0; $j < @ID_to_method; $j++){
	$val = 0;
	push(@tab, $val);
    }
    push(@matrix, \@tab);
}


#Write the pairwise comparison matrises
my $inter;
my $matrix_file;
my %method_sample_driver = ();
foreach $threshold_type (keys %threshold_list){
    foreach $threshold_val (@{$threshold_list{$threshold_type}}){
	clear_matrix(\@matrix);
	foreach $method_1 (@ID_to_method){
	    #Init the mehtod_driver_sample structure
	    my %map = ();
	    $method_sample_driver{$method_1} = \%map;
	    #
	    foreach $method_2 (@ID_to_method){
		#next if($method_1 eq $method_2);
		$inter = 0;
		$list_size = 0;
		#print STDERR $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
		foreach $gene (keys %{$method_result{$method_1}}){
		    
		    #To test the gene status selection
		    next if(! exists $data_gene_annotation{$gene} || #should be replaced
			    ($gene_status_selection eq "CANCER" && $data_gene_annotation{$gene}->{"STATUS"} ne $gene_status_selection));

		    #The must pass the threshold for method_1
		    #print $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
		    if(test_threshold($threshold_type, $threshold_val, $gene, $method_1)){
			
			#Single method analysis
			if($method_1 eq $method_2){
			    $sample_list = $data_gene_annotation{$gene}->{"SAMPLE"};
			    foreach $sample (@{$sample_list}){
				if(! exists $method_sample_driver{$method_1}->{$sample}){
				    $method_sample_driver{$method_1}->{$sample} = 0;
				}
				$method_sample_driver{$method_1}->{$sample}++;
			    }
			}
			#Pairwise comparisions
			else{
			    $list_size++;
			    $method_result{$method_1}->{$gene}->{"USE"} = 1;
			    
			    #print STDERR $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
			    if(exists $method_result{$method_2}->{$gene} && test_threshold($threshold_type, $threshold_val, $gene, $method_2)){
				$inter++;
				$method_result{$method_1}->{$gene}->{"FREQ_CALL"}++;
			    }
			}
		    }
		}
		$matrix[$method_ID{$method_1}]->[$method_ID{$method_2}] = sprintf("%.2f", $inter/$list_size) if($list_size != 0);
	    }
	}
	
	#For sample driver coverage boxplot files
	$boxplot_file = write_boxplot_file($threshold_type, $threshold_val, $gene_status_selection);
	plot_boxplot($boxplot_file, "$threshold_type\_$threshold_val");

	#For the barplot pairwise comparision files
	#$matrix_file = write_barplot_file($threshold_type, $threshold_val, $gene_status_selection);
	#plot_barplot($matrix_file, "$threshold_type\_$threshold_val", @ID_to_method+0, 1);
	#plot_barplot($matrix_file, "$threshold_type\_$threshold_val", @ID_to_method+0, 0);
	
	####For the HEAT map pairwise comparision files
	#Write the file
	#$matrix_file = write_heat_map_file(\@matrix, $threshold_type, $threshold_val, $gene_status_selection);
	#plot_heat_map($matrix_file, "$threshold_type\_$threshold_val");

    }
}

sub write_boxplot_file{
     my ($threshold_type, $threshold_val, $gene_status_selection) = @_;
     my $boxplot_file = "$out_dir/sample_cov_$threshold_type\_$threshold_val";
     if($gene_status_selection ne ""){
	$boxplot_file .= "\_$gene_status_selection";
    }
     
     open(OUT, ">$boxplot_file.dat");

     #The header
     print OUT "".(join("\t", @ID_to_method))."\n";
     
     foreach $sample (keys %sample_list){
	 print OUT $sample;
	 for(my $i = 0; $i < @ID_to_method; $i++){
	     $method = $ID_to_method[$i];
	     $res = 0;
	     if(exists $method_sample_driver{$method}->{$sample}){
		 $res = $method_sample_driver{$method}->{$sample};
	     }
	     print OUT "\t".$res;
	 }
	 print OUT "\n";
     }
     
     close(OUT);

     return $boxplot_file;

}

sub plot_boxplot{
    my ($matrix_file, $title) = @_;
    
    my $font_size = 3;
    my $note_font_size = 8;
    my $margin_size = 30; 

    open(OUT, ">$matrix_file.R");
    print OUT "pdf(file=\"$matrix_file.pdf\",
	paper=\"special\",
	width=10,
	height=10
	)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";

    print OUT "palette <- rainbow(ncol(profile))\n";
    print OUT "boxplot.matrix(as.matrix(profile), col=palette, cex.axis=$font_size)\n";

    close(OUT);
    run_exe("R --vanilla < $matrix_file.R");
}


#Plot the barplot
sub write_barplot_file{
    my ($threshold_type, $threshold_val, $gene_status_selection) = @_;

    #Contruct the frequency call matrix
    #Init the call frequence matrice method
    my @matrix_freq_call = ();
    for(my $i = 0; $i < @ID_to_method; $i++){
	#Init the value 
	$method = $ID_to_method[$i];
	my @tab = ();
	for(my $j = 0; $j < @ID_to_method; $j++){
	    $val = 0;
	    push(@tab, $val);
	}

	foreach $gene (keys %{$method_result{$method}}){
	    #Update the matrix
	    if($method_result{$method}->{$gene}->{"USE"}){
		$tab[$method_result{$method}->{$gene}->{"FREQ_CALL"}]++;
		#Delete the gene value for the next test
		$method_result{$method}->{$gene}->{"FREQ_CALL"} = 0;
		$method_result{$method}->{$gene}->{"USE"} = 0;
	    }
	}

	push(@matrix_freq_call, \@tab);
    }

    my $matrix_file = "$out_dir/freq_call_$threshold_type\_$threshold_val";
    if($gene_status_selection ne ""){
	$matrix_file .= "\_$gene_status_selection";
    }
    open(OUT, ">$matrix_file.dat");
    #header
    print OUT "".join("\t", @ID_to_method)."\n";
    for(my $j = 0; $j < @{$matrix_freq_call[0]}; $j++){
	print OUT $j;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    print OUT "\t".$matrix_freq_call[$i]->[$j];
	}
	print OUT "\n";
    }
    close(OUT);

    return $matrix_file;

}

sub plot_barplot{
    my ($matrix_file, $title, $nb_method, $raw_val) = @_;
    
    my $font_size = 3;
    my $note_font_size = 8;
    my $margin_size = 30; 


    $out_file = $matrix_file;
    $out_file .= "_RAW" if($raw_val);

    open(OUT, ">$out_file.R");
    
    print OUT "pdf(file=\"$out_file.pdf\",
	paper=\"special\",
	width=10,
	height=10
	)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";

    #Compute the freqency
    if(! $raw_val){
	foreach $method (@ID_to_method){
	    print OUT "profile\$$method = profile\$$method/sum(profile\$$method)\n";
	}
    }

    print OUT "profile_mat = as.matrix(profile)\n";
    print OUT "palette <- colorRampPalette(c('#f0f3ff','#0033BB'))($nb_method)\n";
    print OUT "barplot(profile_mat, main=\"$title\", col=palette,  cex.lab=$font_size, cex.names=2, cex.axis = $font_size, cex.main= $font_size)\n";
    
    if(index($matrix_file, "RANK") == -1){
	print OUT "legend(\"top\", legend = rownames(profile), fill = palette, cex=2);"
    }

    close(OUT);
    
    run_exe("R --vanilla < $out_file.R");
}

#Plot the heatmap
sub write_heat_map_file{
    my ($mat, $threshold_type, $threshold_val, $gene_status_selection) = @_;
    my $matrix_file = "$out_dir/pairwise_comparision_$threshold_type\_$threshold_val";
    if($gene_status_selection ne ""){
	$matrix_file .= "\_$gene_status_selection";
    }
    open(OUT, ">$matrix_file.dat");
    #header
    print OUT "".join("\t", @ID_to_method)."\n";
    for(my $i = 0; $i < @ID_to_method; $i++){
	print OUT $ID_to_method[$i];
	for(my $j = 0; $j < @ID_to_method; $j++){
	    print OUT "\t".$mat->[$i]->[$j];
	}
	print OUT "\n";
    }
    close(OUT);
    
    return $matrix_file;

}

sub plot_heat_map{
    my ($matrix_file, $title) = @_;

    my $font_size = 6;
    my $note_font_size = 8;
    my $margin_size = 30; 

    open(OUT, ">$matrix_file.R");
    
    print OUT "pdf(file=\"$matrix_file.pdf\",
	paper=\"special\",
	width=15,
	height=15
	)\n";

    print OUT "library(\"gplots\")\n";
    print OUT "palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(10)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";
    print OUT "profile_mat = as.matrix(profile)\n";
#print OUT "heatmap.2(profile_mat, col=redgreen, margin=c(5, 20), key=TRUE, scale=\"row\", density.info=\"none\", trace=\"none\")\n";
    
    
    my $str_key = "key = FALSE,
lmat=rbind(c(2),c(3),c(1),c(3)), 
    lhei=c(1,1,9,0), 
    lwid=c(1)";
    
    if($matrix_file eq "$out_dir/pairwise_comparision_MUT_FREQ_0" || $matrix_file eq "$out_dir/pairwise_comparision_MUT_FREQ_0_CANCER"){
	$str_key = "key=TRUE, keysize=1.3";
	$note_font_size = 4;
	$font_size = 6;
	$margin_size = 30; 
    }
    
    print OUT 
	"heatmap.2(profile_mat,
main = \"$title\", 
scale=\"none\",
Rowv=FALSE,Colv=FALSE,
breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), #symbreaks = TRUE,
#
#key.par=list(mgp=c(1.5, 0.5, 0),
#mar=c(2.5, 2.5, 1, 0)),
#
dendrogram = \"none\", 
$str_key,
cellnote=as.matrix(profile_mat),notecol=\"black\",notecex=$note_font_size,
               #hclustfun = function(x) hclust(x,method = 'complete'),
               #distfun = function(x) dist(x,method = 'euclidean'),
               margin=c($margin_size, $margin_size), 
col=palette, cexRow=$font_size, cexCol=$font_size, 
               density.info=\"none\", trace=\"none\"";
    #print OUT ",ColSideColors = $color_subtype" if($subtype_file ne "NONE");
    print OUT ")\n";
    close(OUT);

    run_exe("R --vanilla < $matrix_file.R");
    run_exe("pdfcrop  $matrix_file.pdf $matrix_file\_temp.pdf");
    run_exe("mv  $matrix_file\_temp.pdf $matrix_file.pdf");

}

sub cancer_annotation{
    #cancer gene census
    my $cancer_gene_census_file = "/mnt/pnsg10_projects/bertrandd/oncoimpact/SCRIPT/oncoIMPACT/cancer_gene_census.csv";
    open(FILE, $cancer_gene_census_file);
    while(<FILE>){
	@line = split(/\s+/, $_);    
	#chomp($line[0]);
	$cancer_census{$line[0]} = 1;
    }
    close(FILE);

#pan cancer data
    my $pan_cancer_file = "/mnt/pnsg10_projects/bertrandd/oncoimpact/SCRIPT/oncoIMPACT/pancancer_driver_list.csv";
    open(FILE, $pan_cancer_file);
    while(<FILE>){
	chop($_);
	@line = split(/\t/, $_);
	#chomp($line[0]);
	$pan_cancer{$line[0]} = 1 if($line[7] eq "High_Confidence_Driver");
    }
    close(FILE);
}


sub clear_matrix{
    my ($mat) = @_;
    for(my $i = 0; $i < @{$mat}; $i++){
	for(my $j = 0; $j < @{$mat->[$i]}; $j++){
	    $val = 0;
	    $val = 1 if($i == $j); 
	    $mat->[$i]->[$j] = $val;
	}
    }
}

sub test_threshold{
    my ($t_type, $t_val, $gene, $meth) = @_;
    if(#Gene annotation base threshold
       ($t_type =~ m/MUT_FREQ/ && $data_gene_annotation{$gene}->{$t_type} >= $t_val) ||
       #Method output based threshold
       ($t_type =~ m/RANK/ && $method_result{$meth}->{$gene}->{$t_type} <= $t_val)){
	return 1;
    }
    return 0;
}


sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";
    print STDERR `$exe` if($run);
}
