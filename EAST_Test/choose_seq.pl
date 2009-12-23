#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;

my $min_length = 0;
my $max_length = 100000;
my $index = -1;

my $marker = 0;
while (substr($ARGV[$marker], 0, 1) eq '-') {
    $ARGV[$marker] =~ s/^\-+//;
    $min_length = $ARGV[++$marker] if (grep {$ARGV[$marker] eq $_} qw(min min_length));
    $max_length = $ARGV[++$marker] if (grep {$ARGV[$marker] eq $_} qw(max max_length));
    $index = $ARGV[++$marker] if (grep {$ARGV[$marker] eq $_} qw(index i));
    $marker++;
}

my ($input, $output) = @ARGV[$marker, $marker+1];

my $seqin = Bio::SeqIO->new(-format => "fasta", -file => $input);
my $seqout = Bio::SeqIO->new(-format => "fasta", -file => ">".$output);

my @arr;
while (my $seq = $seqin->next_seq) {
    push @arr, $seq if $seq->length >= $min_length && $seq->length <= $max_length
}

$index = int(rand(@arr)) if $index == -1;
print scalar @arr, "\n", $index, "\n";
$seqout->write_seq($arr[$index]);
