#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <trinity_cdhit.fa> <transdecoder.fa>\n" unless @ARGV == 2;

my $trinity = shift;
my $transdecoder = shift;


my %hash;
open IN, $trinity || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$hash{$1} = 1;
	}
}
close IN;

my (%seq, $name);
open TT, $transdecoder || die "$!\n";
while( <TT> ){
	chomp;
	if( />(\S+)/ ){
		$name = $1;
	}else{
		$seq{$name} .= $_;
	}
}
close TT;

for my $key ( sort keys %seq ){
	my $tt = (split(/\./, $key))[0];
	if( exists $hash{$tt} ){
		print ">$key\n$seq{$key}\n";
	}
}
