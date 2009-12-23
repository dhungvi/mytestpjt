#!/usr/bin/perl -w
use strict;
use POSIX qw(log10);

use Bio::SeqIO;
my $filename  = $ARGV[0];
my $errorRate = $ARGV[1];

if ( $errorRate != 0 ) {
	#my $quality = int( -10 * ( log($errorRate) / log(10) ) );
	my $quality = int(-10 * (log10($errorRate)));
	my $seq_in = Bio::SeqIO->new( -format => "fasta", -file => "$filename" );
	open OFILE, ">$filename.qual"
	  or die("Fail to open the output file $filename.qual!\n");

	while ( my $seq_obj = $seq_in->next_seq ) {
		my $id = $seq_obj->id;
		print OFILE ">$id\n";
		my $len = length( $seq_obj->seq );
		my $qual;
		my $num = 0;
    	for (my $i=0; $i < $len; $i++) {
    		if ($num == 20) {
    			$qual .= "\n";
    			$num = 0;
    		}
 			$qual .= "$quality ";
    		$num++;
		}

		print OFILE "$qual\n";
	}

	close OFILE;
}
