#!/usr/bin/perl -w
use strict;
use redhawk;
use Bio::SeqIO;


# Parameters relationg to redhawk
my $job_limit = 500;
my $pause = 10;

# Parameters relating to actual testing
my $geneFile;
my $geneSource = "all_zf_cdnas.reduced.fa";
my $num_genes = $ARGV[0];
my $minGeneLen = $ARGV[1];
my $maxGeneLen = $ARGV[2];
my $useQualityFile = $ARGV[3]; #1-use; 0-not use.
my $numOfTrials = 60;

my @gene_list = choose_genes($geneSource, $num_genes, $minGeneLen, $maxGeneLen);

# In the following we execute all batch files, and keep a list in the @batch_files array
my @batch_files;
foreach my $geneSeq (@gene_list) {
    my $str = $geneSeq->seq;
    my $len = length($str);
    my $numOfEst = int(25*$len/360);

    # Write gene to a file
    my $geneID = $geneSeq->id;
    my $geneFile = "geneFile.$geneID.$$.fa";
    my $seqout = Bio::SeqIO->new(-format => "fasta", -file => ">$geneFile");
    $seqout->write_seq($geneSeq);


    for ( my $errorRate = 0 ; $errorRate <= 8; $errorRate += 0.25) {
	my $rate = $errorRate / 100;
	
	for ( my $i = 0 ; $i < $numOfTrials ; $i++ ) {
	    my $uniqueId = "$$\_$geneID\_$errorRate\_$i";
	    my $executable = "module load cap3; module load bioperl; ./singleTest.pl $geneFile $numOfEst $rate $uniqueId $useQualityFile"; 
	    
	    # Set up readhawk execution       
	    my $file = create_batch_file(file => "singleTest.$uniqueId", executable => $executable);
	    my $id = execute_batch($file);
	    push(@batch_files, [$file,$id,$i,$errorRate,$geneID,$len,$numOfEst]);
	    wait_on_job_limit($job_limit, $pause);
	}
    }

    #system("rm $geneFile");   # Clean up the files (do we want them cleaned up??? no!)
}

# Now we go through the @batch_files array, process the results and write them to a file.

#Open file and print comments (as before)
use POSIX qw( strftime );
my $yyyymmdd = strftime( "%Y-%m-%d %H-%M-%S", localtime );    #get current time

my $output_file = "results\_$num_genes\_$minGeneLen\_$maxGeneLen\_$yyyymmdd";

open OFILE, ">$output_file.txt" or die ("Fail to open the output file $output_file.txt!\n"); 
open EFILE, ">$output_file.err" or die("Fail to open the eroor file $output_file.err\n");
my $output = join "\t", qw(geneID geneLength numberOfEsts ErrorRate Cap3NumContigs EastNumContigs Cap3NumSing EastNumSin Cap3Len EastLen Cap3AScore EastAScore Cap3Time EastTime Cap3NumUsedEsts EastNumUsedEsts);
print OFILE "#number_of_runs=$numOfTrials\tcoverage=25\tuse_quality_file=$useQualityFile\n#Date:$yyyymmdd\n$output\n";

foreach my $p (@batch_files) {
    my ($file, $id, $trial, $errorRate, $geneID, $geneLen, $numEST) = @$p;
    my $rate = $errorRate / 100;

    wait_on_job($file, $id);

    # Read in output messages
    my $fh = ofile_handle($file, $id);
    my $output = <$fh>;
    chomp $output;
    my @results = split(/\s+/, $output);
    if (@results == 12) {
	print OFILE join "\t", ($geneID, $geneLen, $numEST, $rate);
	print OFILE "\t", $output;
	print OFILE "\n";
    }
    close($fh);
    
    # Read in error messages
    my $efh = efile_handle($file, $id);
    my @input = <$efh>;
    if (@input) {
	print EFILE "gene: $geneID, Error: $errorRate, trial: $trial\n";
	print EFILE join "\n", @input;
	print EFILE "\n";
    }
    
    cleanup($file,$id);
}

sub choose_genes {
    my $file = $_[0];
    my $num_seq = $_[1];
    my $min_length = $_[2] || 0;
    my $max_length = $_[3] || 10000000;

    my @results;

    my $seqin = Bio::SeqIO->new(-format => "fasta", -file => $file);
    
    my @arr;
    while (my $seq = $seqin->next_seq) {
	push @arr, $seq if $seq->length >= $min_length && $seq->length <= $max_length
    }
    
    return @arr if $num_seq == -1;
    die("choose_seq: Num enough sequence in file\n") if (@arr < $num_seq);
    for (my $i=0; $i < $num_seq; $i++) {
	my $index = int(rand(@arr));
	push @results, $arr[$index];
	$arr[$index] = pop(@arr);
    }

    return @results;
}

chdir "/home/zhangy9/peace";
system("rm *.job");
system("rm *.o*");
