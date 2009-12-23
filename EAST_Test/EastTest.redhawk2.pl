#!/usr/bin/perl -w
use strict;
use redhawk;
use Bio::SeqIO;

# Parameters relationg to redhawk
my $job_limit = 100;
my $pause = 10;

# Parameters relating to actual testing
my $geneFile;
my $minGeneLen = $ARGV[0];
my $maxGeneLen = $ARGV[1];
my $numOfTrials = 60;

# Randomly take a gene
if (! defined $ARGV[1]) {
	$geneFile = "geneFile\_$minGeneLen.fa";
	system("perl choose_seq.pl -min $minGeneLen all_zf_cdnas.fa $geneFile");
} else {
	$geneFile = "geneFile\_$minGeneLen\_$maxGeneLen.fa";
	system("perl choose_seq.pl -min $minGeneLen -max $maxGeneLen all_zf_cdnas.fa $geneFile");
}
system("cp -f $geneFile ~/peace/");
my $seqin = Bio::SeqIO->new(-format => 'fasta', -file => $geneFile);
my $str = $seqin->next_seq->seq;
my $len = length($str);
my $numOfEst = int(25*$len/360);

my %batch_files;
for ( my $error_rate = 0 ; $error_rate <= 8 ; $error_rate += 0.25 ) {
    $batch_files{$error_rate} = [];
    my $rate = $error_rate / 100;
    
    for ( my $i = 0 ; $i < $numOfTrials ; $i++ ) {
	my $uniqueId = "$error_rate\_$i";
	my $executable = "module load cap3; module load bioperl; ./singleTest.pl $geneFile $numOfEst $rate $uniqueId 0"; 

	# Set up readhawk execution       
	my $file = create_batch_file(file => "singleTest.$uniqueId", executable => $executable);
	my $id = execute_batch($file);
	push(@{$batch_files{$error_rate}}, [$file,$id, $i]);
	wait_on_job_limit($job_limit, $pause);
    }
}


#Open file and print comments (as before)
use POSIX qw( strftime );
my $yyyymmdd = strftime( "%Y-%m-%d %H-%M-%S", localtime );    #get current time

my $output_file = "$geneFile\_$yyyymmdd";

open OFILE, ">$output_file.txt" or die ("Fail to open the output file $output_file.txt!\n"); 
open EFILE, ">$output_file.err" or die("Fail to open the eroor file $output_file.err\n");
my $output = "Cap3NumContigs\tEastNumContigs\tCap3NumSing\tEastNumSin\tCap3Len\tEastLen\tCap3AScore\tEastAScore\tCap3Time\tEastTime\tCap3NumUsedEsts\tEastNumUsedEsts\tErrorRate";
print OFILE "#$geneFile, length_of_gene=$len, number_of_ests=$numOfEst, number_of_runs=$numOfTrials\n#Date:$yyyymmdd\n$output\n";

for my $error_rate (sort {$a <=> $b} keys %batch_files) {
    my $rate = $error_rate / 100;

    #Extracts information output information
    for my $pair (@{$batch_files{$error_rate}}) {
	my ($file, $id, $trial) = @$pair;
	wait_on_job($file, $id);

	# Read in output messages
	my $fh = ofile_handle($file, $id);
	my $output = <$fh>;
	chomp $output;
	my @results = split(/\s+/, $output);
	print OFILE "$output\t$rate\n" if @results == 12;
	close($fh);

	# Read in error messages
	my $efh = efile_handle($file, $id);
	my @input = <$efh>;
	if (@input) {
	    print EFILE "Error: $error_rate, trial: $trial\n";
	    print EFILE join "\n", @input;
	    print EFILE "\n";
	}

	cleanup($file,$id);
    }
}

