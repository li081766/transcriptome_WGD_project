#!/usr/bin/perl -w
use strict;
die "perl $0 <in.fasta> <cutoff>\n\tdefault cutoff is 150bp\n" if @ARGV < 1;

my $in = shift;

my $cutoff;
if( @ARGV == 2 ){
	$cutoff = shift;
}else{
	$cutoff = 150;
}

open IN, $in || die "$!\n";
my ($id, %seq);
while( <IN> ){
	chomp;
	if( />(\S+).*/ ){
		$id = $1;
	}else{
		$seq{$id} .= $_
	}
}
close IN;

for my $i ( sort keys %seq ){
	if( length($seq{$i}) >= $cutoff ){	
		print ">$i\n$seq{$i}\n";
	}
}
