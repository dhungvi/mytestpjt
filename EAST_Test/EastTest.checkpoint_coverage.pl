#!/usr/bin/perl -w
use strict;
use redhawk;
use Bio::SeqIO;

# Parameters relationg to redhawk
my $job_limit = 100;
my $pause     = 10;

# Parameters relating to actual testing
my $job_type = 0;  # 1 = execute only; 2 = process only; 0 = both
my $reset = 0;
while ($ARGV[0] =~ /^-/) {
    my $switch = shift @ARGV;
    $switch =~ s/^-+//;
    if (grep {$switch eq $_} qw(e execute)) {
	$job_type = 1;
    }
    if (grep {$switch eq $_} qw(p process)) {
	$job_type = 2;
    }
    if (grep {$switch eq $_} qw(b both)) {
	$job_type = 0;
    }
    if (grep {$switch eq $_} qw(jt job_type)) {
	$job_type = shift @ARGV;
    }
    if (grep {$switch eq $_} qw(r reset)) {
	$reset = 1;
    }
}

my $batch_record_file = $ARGV[0] || die("Need a batch_record_file");
my $progress_file = $ARGV[1];

my $geneFile;
my $geneSource     = "all_zf_cdnas.reduced.fa";
my $num_genes      = $ARGV[2];
my $minGeneLen     = $ARGV[3];
my $maxGeneLen     = $ARGV[4];
my $useQualityFile = $ARGV[5];                   #1-use; 0-not use.
my $numOfTrials    = 60;

my @gene_list;
my @batch_files;
my %processed_jobs;
# Read in record
if (-e $batch_record_file && $reset == 0) {
    die("Batch record but no gene record\n") unless -e "$batch_record_file.fa";
    my $seqin = Bio::SeqIO->new(-format => "fasta", -file => "$batch_record_file.fa");
    while (my $seq = $seqin->next_seq) {
	push(@gene_list, $seq);
    }

    open(BFILE, "$batch_record_file");
    while (<BFILE>) {
	chomp;
	my @batch_info = split "\t", $_;
	push @batch_files, \@batch_info;
	my $key = join "\t", @batch_info[2..(@batch_info-1)];
	$processed_jobs{$key} = 1;
    }
    close(BFILE);
}
else {
    @gene_list = choose_genes( $geneSource, $num_genes, $minGeneLen, $maxGeneLen );
	my $seqout = Bio::SeqIO->new(-format => "fasta", -file => ">$batch_record_file.fa");
    foreach(@gene_list) {
	$seqout->write_seq($_);
    }
    unlink "$batch_record_file";
}

open(PFILE, ">$progress_file");

if ($job_type != 2) {
    open(BFILE, ">>$batch_record_file");

# In the following we execute all batch files, and keep a list in the @batch_files array
my $gene_count = 0;
foreach my $geneSeq (@gene_list) {
	my $str = $geneSeq->seq;
	my $len = length($str);
	my $numOfEst;

	# Write gene to a file
	my $geneID   = $geneSeq->id;
	my $geneFile = "geneFile.$geneID.$$.fa";
	my $seqout   = Bio::SeqIO->new( -format => "fasta", -file => ">$geneFile" );
	$seqout->write_seq($geneSeq);

	for ( my $errorRate = 1 ; $errorRate <= 6 ; $errorRate += 2.5 )
	{
		for ( my $coverage = 5 ; $coverage <= 50 ; $coverage += 5 ) {
			print PFILE "Starting jobs: gene $gene_count, error rate $errorRate, coverage $coverage\n";
			my $rate = $errorRate / 100;
			$numOfEst = int( $coverage * $len / 360 );

			for ( my $i = 0 ; $i < $numOfTrials ; $i++ ) {
			my @batch_info = ($i, $rate, $coverage, $geneID, $len, $numOfEst);
			unless ($processed_jobs{join "\t", @batch_info}) {

				my $uniqueId   = "$$\_$geneID\_$errorRate\_$coverage\_$i";
				my $executable = "module load cap3; module load bioperl; ./singleTest.pl $geneFile $numOfEst $rate $uniqueId $useQualityFile";

				# Set up readhawk execution
				my $file = create_batch_file(
					file       => "singleTest.$uniqueId",
					executable => $executable
				);
				my $id = execute_batch($file);
				push(@batch_files,[$file, $id, @batch_info]);
				print BFILE join "\t", @{$batch_files[-1]}, "\n";
				
				wait_on_job_limit( $job_limit, $pause );
			}
			}
		}
		$gene_count = $gene_count + 1;
	}

#system("rm $geneFile");   # Clean up the files (do we want them cleaned up??? no!)
}
close(BFILE);
}

# Now we go through the @batch_files array, process the results and write them to a file.
if ($job_type != 1) {

#Open file and print comments (as before)
use POSIX qw( strftime );
my $yyyymmdd = strftime( "%Y-%m-%d %H-%M-%S", localtime );    #get current time

my $output_file =
  "resultsOnCoverage\_$num_genes\_$minGeneLen\_$maxGeneLen\_$yyyymmdd";

open OFILE, ">$output_file.txt"
  or die("Fail to open the output file $output_file.txt!\n");
open EFILE, ">$output_file.err"
  or die("Fail to open the eroor file $output_file.err\n");
my $output = join "\t",
  qw(geneID geneLength numberOfEsts ErrorRate Coverage Cap3NumContigs EastNumContigs Cap3NumSing EastNumSin Cap3Len EastLen Cap3AScore EastAScore Cap3Time EastTime Cap3NumUsedEsts EastNumUsedEsts);
print OFILE
"#number_of_runs=$numOfTrials\tuse_quality_file=$useQualityFile\n#Date:$yyyymmdd\n$output\n";

foreach my $p (@batch_files) {
	my ( $file, $id, $trial, $rate, $coverage, $geneID, $geneLen, $numEST ) = @$p;

	wait_on_job( $file, $id );

	# Read in output messages
	my $fh = ofile_handle( $file, $id );
	my $output = <$fh>;
	chomp $output;
	my @results = split( /\s+/, $output );
	if ( @results == 12 ) {
		print OFILE join "\t", ( $geneID, $geneLen, $numEST, $rate, $coverage );
		print OFILE "\t", $output;
		print OFILE "\n";
	}
	close($fh);

	# Read in error messages
	my $efh = efile_handle( $file, $id );
	my @input = <$efh>;
	if (@input) {
		print EFILE "gene: $geneID, errorRate: $rate, coverage: $coverage, trial: $trial\n";
		print EFILE join "\n", @input;
		print EFILE "\n";
	}

	cleanup( $file, $id );
}
}

sub choose_genes {
	my $file       = $_[0];
	my $num_seq    = $_[1];
	my $min_length = $_[2] || 0;
	my $max_length = $_[3] || 10000000;

	my @results;

	my $seqin = Bio::SeqIO->new( -format => "fasta", -file => $file );

	my @arr;
	while ( my $seq = $seqin->next_seq ) {
		push @arr, $seq
		  if $seq->length >= $min_length && $seq->length <= $max_length;
	}

	return @arr                                      if $num_seq == -1;
	die("choose_seq: Num enough sequence in file\n") if ( @arr < $num_seq );
	for ( my $i = 0 ; $i < $num_seq ; $i++ ) {
		my $index = int( rand(@arr) );
		push @results, $arr[$index];
		$arr[$index] = pop(@arr);
	}

	return @results;
}

chdir "../peace";
system("rm *.job");
system("rm *.o*");
