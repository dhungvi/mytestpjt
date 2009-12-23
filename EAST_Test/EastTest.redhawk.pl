#!/usr/bin/perl -w
use strict;
use redhawk;

# Parameters relationg to redhawk
my $job_limit = 100;
my $pause = 10;

# Parameters relating to actual testing
my $geneFile    = "geneFile.fa";
my $numOfEst    = 300;
my $numOfTrials = 10;



my %batch_files;
for ( my $error_rate = 0 ; $error_rate <= 1 ; $error_rate += 0.25 ) {
    $batch_files{$error_rate} = [];
    my $rate = $error_rate / 100;
    
    for ( my $i = 0 ; $i < $numOfTrials ; $i++ ) {
	my $uniqueId = "$error_rate\_$i";
	my $executable = "module load cap3; module load bioperl; ./singleTest.pl $geneFile $numOfEst $rate $uniqueId"; 

	# Set up readhawk execution       
	my $file = create_batch_file(file => "singleTest.$uniqueId", executable => $executable);
	my $id = execute_batch($file);
	push(@{$batch_files{$error_rate}}, [$file,$id]);
	wait_on_job_limit($job_limit, $pause);
    }
}


for my $error_rate (sort {$a <=> $b} keys %batch_files) {
    my $rate = $error_rate / 100;

    #Open file and print comments (as before)
    open OFILE, ">test.$error_rate.txt" or die ("Fail to open the output file $error_rate.txt!"); 
    my $output = "#Cap3NumContigs  EastNumContigs  Cap3NumSing  EastNumSin  Cap3Len  EastLen  Cap3AScore  EastAScore  Cap3Time  EastTime  Cap3NumUsedEsts  EastNumUsedEsts";
    use POSIX qw( strftime );
    my $mmddyyyy = strftime( "%m/%d/%Y", localtime );    #get current time
    print OFILE "$output\n#error_rate=$rate, number_of_ests=$numOfEst, number_of_runs=$numOfTrials\n#Date:$mmddyyyy\n";

    #Extracts information output information
    for my $pair (@{$batch_files{$error_rate}}) {
	my ($file, $id) = @$pair;
	wait_on_job($file, $id);
	my $fh = ofile_handle($file, $id);
	my $output = <$fh>;
	chomp $output;
	print OFILE "$output\n";
	cleanup($file,$id);
    }
    close OFILE;
}

