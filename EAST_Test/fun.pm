#!/usr/bin/perl -w
use strict;

package fun;

######################
# Author: John Karro
# Date: 2/6/06
# 
# Useful Perl functions.
######################

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(next_line ckopen copen convert_nucleotide min max fopen vectorPair2Hash org_to_build build_to_org org_and_build is_org is_build build_to_chr_array stats parse_location parse_cmdln);
our $VERSION = 1.00;

#############################################
# next_line:
#
# Return the next non-blank line in the file
#############################################
sub next_line {
    my $fh = $_[0];
    $fh || die("No file pointer provided to next_line\n");

    my $line;
    while ($line = <$fh>) {
	chomp($line);
	return $line unless $line =~ /^\s*$/;
    }

    return undef;
}
	

#######################################################
# convert_nucleotide:
# 
# Convert a nucleotide into number (0-4, alphabetical)
#######################################################
sub convert_nucleotide {
    my $n = $_[0];

    return -1 if $n eq '-';
    return  0 if $n eq 'a' || $n eq 'A';
    return  1 if $n eq 'c' || $n eq 'C';
    return  2 if $n eq 'g' || $n eq 'G';
    return  3 if $n eq 't' || $n eq 'T';
    return -2;
}


###########################################
# copen
#
# Open a file and return the file handle.
# Dies and prints file name of file cannot
# be opened.
###########################################
sub ckopen {
    my $file = $_[0];
    $file || die("copen requires a file name\n");

    return *STDIN if $file eq '-';

    my $fp;
    open($fp, $file) || die("Can't open file: $file\n");
    
    return $fp;
}

sub copen {
    return ckopen @_;
}

#########################################
# fopen
#
# Open a file and return the file handle.
# Returns undef if file cannot be opened.
#########################################
sub fopen {
    my $file = $_[0];

    return *STDIN if $file eq '-';

    my $fp;
    open($fp, $file) || return undef;

    return $fp;
}


#########################################
# min / max
#
# Return the min / max of two numbers
#########################################
sub min {
    return $_[0] <= $_[1] ? $_[0] : $_[1];
}

sub max {
    return $_[0] > $_[1] ? $_[0] : $_[1];
}

#########################################
# vectorPair2Hash
# 
# Take two vector refs of the same length and
# return a hash matching corresponding
# elements, first to second.
#########################################
sub vectorPair2Hash {
    my @arr1 = @{$_[0]};
    my @arr2 = @{$_[1]};
    
    @arr1 == @arr2 || die("vectorPair2Hash: Vectors of differing lengths\n");

    my %H;
    @H{@arr1} = @arr2;
  
    return %H;
}

#########################################
# trim
# 
# Remove end spaces
#########################################
sub trim {
    $_[0] =~ s/^\s+//;
    $_[0] =~ s/\s+$//;
    return $_[0];
}

#########################################
# org_to_build
#
# Return the current build name for a 
# given organism (or undef)
#########################################
my %org_to_build = ("human" => "hg18", 
		    "chimp" => "pt2", 
		    "macaque" => "rm2",
		    "mouse" => "mm8",
		    "rat" => "rn3",
		    "dog" => "cf2",
		    "chicken" => "gg2",
		    "fly" => "dm3",
		    "worm" => "ce3");

my %build_to_org = reverse(%org_to_build);
sub org_to_build {
    return $org_to_build{$_[0]} || undef;
}

#########################################
# build_to_org
#
# Returns the organism name for a known
# build (or undef)
#########################################
sub build_to_org {
    return $build_to_org{$_[0]} || undef;
}

#########################################
# build_and_org
#
# Take a build or org and return
# the build and org pair
#########################################
sub org_and_build {
    my $build = org_to_build($_[0]) || $_[0];
    my $org = build_to_org($build);

    return $org ? ($org, $build) : (undef, undef);
}


#########################################
# is_org
#
# Returns true if parameter is a known
# org. name.
#########################################
sub is_org {
    return exists $org_to_build{$_[0]};
}

#########################################
# is_build
#
# Returns true if parameter is a known
# build.
#########################################
sub is_build {
    return exists $build_to_org{$_[0]};
}

#########################################
# build_to_chr_array
#
# Return th chromosome list for a given
# build.  Second parameter indicates
# a request for autosomes only.
#########################################
my %build_to_chr_array = ('hg18' => [[1..22], ['X', 'Y']],
		'cf2' => [[1..38], ['X']],
		'mm8' => [[1..19], ['X',  'Y']],
		'rn3' => [[1..20], ['X']],
		'pt2' => [[1, '2a', '2b', 3..22], ['X', 'Y']],
		'rm2' => [[1..20], ['X']],
		'gg2' => [[1..24, 26..28, 32], ['W', 'Z']],
		'canHg12' => [[1..30]],
                'dm3' => [['2L', '2LHet', '2R', '2RHet', '3L', '3LHet', '3R', '3RHet', '4']],
                'ce4' => [['I', 'II', 'III', 'IV', 'V']]);

sub build_to_chr_array {
    my $r = $build_to_chr_array{$_[0]};
    my $auto_only = $_[1];
    $r || die("fun::build_to_chr_array: bad org($_[0])\n");

    my @chr = $auto_only ? @{$r->[0]} : (@{$r->[0]}, @{$r->[1]});

    return @chr;
}

########################################
# stats
# Return mean and var. of the arguments
########################################
sub stats {
    my @arr = @_;
    @arr > 0 || die("fun::stats: empty argument\n");
   
    my $total = 0;
    for (my $i=0; $i < @arr; $i++) {
	$total += $arr[$i];
    }
    my $mean = (1.0*$total)/@arr;

    $total = 0;
    for (my $i=0; $i < @arr; $i++) {
	my $d = $arr[$i]-$mean;
	$total += $d*$d;
    }

    my $var = $total / (@arr-1);

    return ($mean, $var);
}

#########################################
# parse_location
#
# Parse a chromsome location of the form
# chrX:a-b (returns (chrX, a, b))
#########################################
sub parse_location {
    return $_[0] =~ /(chr\w+):(\d+)-(\d+)/g;
}


#######################
# parse_cmdln
#
# Parse argv for switches
#######################
sub parse_cmdln {
    my @argv = @{shift @_};
    my %H = @_;

    while (@argv && substr($argv[0], 0, 1) eq '-') {
	my $switch = shift @argv;
	$switch =~ s/^\-\-/\-/;
	my $discriptor_p = $H{$switch} || die("Bad cmd. line switch: $switch\n");
	if (ref($discriptor_p) eq "SCALAR") {
	    $$discriptor_p = shift @argv || die("Missing cmd. line parameter after: $switch\n");
	}
	elsif (ref($discriptor_p) eq "CODE") {
	    my @params = split /\s+/, shift @argv;
	    &{$discriptor_p}(@params);
	}
	elsif (ref($discriptor_p) eq "ARRAY") {	    
	    my @discriptor = @$discriptor_p;

 	    if (ref($discriptor[0]) eq "CODE") {
                my $f = shift @discriptor;
		my @params = map {split /\s+/, $_} @discriptor;
		&{$f}(@params);
	    }
	    else {
		my $p = $discriptor[0];
		if (@discriptor == 1) {
		    $$p = shift @argv || die("Missing cmd. line parameter after: $switch\n");
		}
		else {
		    $$p = $discriptor[1] || die("Bad discriptor given to parse_cmdln\n");
		}
	    }
	}
	else {
	    die("Bad discriptor given to parse_cmdln\n");
	}
    }

    return @argv;
}
