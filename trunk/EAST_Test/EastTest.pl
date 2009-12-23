#!/usr/bin/perl -w
use strict;

my $geneFile    = "geneFile.fa";
my $numOfEst    = 300;
my $numOfTrials = 30;

for ( my $error_rate = 0 ; $error_rate <= 8 ; $error_rate += 0.25 ) {
	my $rate = $error_rate / 100;
	open OFILE, ">test.$error_rate.txt" or die ("Fail to open the output file test.$error_rate.txt!"); 
	my $output = "#Cap3NumContigs  EastNumContigs  Cap3NumSing  EastNumSin  Cap3Len  EastLen  Cap3AScore  EastAScore  Cap3Time  EastTime  Cap3NumUsedEsts  EastNumUsedEsts";
	use POSIX qw( strftime );
	my $mmddyyyy = strftime( "%m/%d/%Y", localtime );    #get current time
	print OFILE "$output\n#error_rate=$rate, number_of_ests=$numOfEst, number_of_runs=$numOfTrials\n#Date:$mmddyyyy\n";

	for ( my $i = 0 ; $i < $numOfTrials ; $i++ ) {
		my $t = $i + 1;
		my $uniqueId = "$error_rate\_$t";

		$output = `./singleTest.pl $geneFile $numOfEst $rate $uniqueId 0`;
		chomp $output;
		print OFILE "$output\n";
	}
	close OFILE;
}

