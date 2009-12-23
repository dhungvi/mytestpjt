#!/usr/bin/perl -w
use strict;

package redhawk;

require Exporter;

# Documentation at end of file

our @ISA = qw(Exporter);
our @EXPORT = qw(create_batch_file execute_batch cleanup ofile_handle efile_handle wait_on_job wait_on_job_limit);
our $VERSION = 1.00;

sub create_batch_file {
    my %H = @_;

    exists $H{file} || die("create_batch_file: Need a file name.\n");
    exists $H{executable} || die("create_batch_file: Need an executable command.\n");
    my $job = $H{name} || $H{file};
    my $nodes = $H{nodes} || 1;
    my $ppn = $H{ppn} || 1;
    my $walltime = $H{walltime} || "10:00:00";

    open(IFILE, ">".$H{file});
    print IFILE "\#\!/bin/bash -l\n";
    print IFILE "\#PBS -N $job\n";
    print IFILE "\#PBS -l nodes=$nodes:ppn=$ppn\n";
    print IFILE "\#PBS -l walltime=$walltime\n";
    print IFILE "\#PBS -M $H{address}\n" if $H{address};
    print IFILE "\#PBS -j $H{join}\n" if $H{join};
    print IFILE "\#PBS -V\n" if $H{env};
    print IFILE "\#PBS -q $H{queue}\n" if $H{queue};
    print IFILE "\#PBS -m $H{mail}\n" if $H{mail};

    
    if ($H{dir}) {
	print IFILE "cd $H{dir}";
    }
    else {
	print IFILE "cd \$PBS_O_WORKDIR\n";
    }
    print IFILE $H{executable}, "\n";

    close(IFILE);
    return $H{file};
}

sub execute_batch {
    my $file = $_[0];
    my $preserve = $_[1];
    my $v = `qsub $file`;
    $v =~ /^(\d+)/;

    unlink $file unless $preserve;
    return $1;
}

sub wait_on_job {
    my ($name, $id) = @_[0,1];
    my $delay = $_[2] || 10;

    my $file = "$name.o$id";
    while (!(-e $file)) {
	sleep($delay);
    }
}

sub ofile_handle {
    my ($name, $id) = @_;
    my $fp;
    open($fp, "$name.o$id");
    return $fp;
}

sub efile_handle {
    my ($name, $id) = @_;
    my $fp;
    open($fp, "$name.e$id");
    return $fp;
}

sub cleanup {
    my ($name, $id) = @_;
    unlink "$name.o$id";
    unlink "$name.e$id";
}

sub wait_on_job_limit {
    my $limit = $_[0];
    my $delay = $_[1] || 10;
    my $user = $_[2] || 'karroje';

    my $d = 0;
    while(1) {
	`qstat -u $user | wc -l` =~ /(\d+)/;
	$d = $1 - 5;
	last if ($d < $limit);
	sleep($delay);
    }
}


