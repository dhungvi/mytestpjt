#!/usr/bin/perl -w
use strict;

my @arr = split(/\n/, `qstat -u zhangy9`);
for (@arr) {
  if (/(\d+)\.mulnx/) {
     my $num = $1;
     `qdel $1`;
  }
}
	
