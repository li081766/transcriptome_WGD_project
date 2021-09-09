#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <cdhit.shell> <out_file>\n" unless @ARGV == 2;

my $in = shift;
my $outfile = shift;

open IN, $in || die "$!\n";
open OUT, ">$outfile\n";
my $num = 0;
while( <IN> ){
	chomp;
	if( /^mv/ ){
		$num++;
		my @tmp = split(/\s+/, $_);
		print OUT ">Cluster $num\n";
		open TT, $tmp[2] || die "$!\n";
		my ($seq, $len);
		while( <TT> ){
			chomp;
			if( /^>/ ){
				$seq = $_;
			}else{
				$len = length($_);
			}
		}
		close TT;
		print OUT "0\t$len", "nt, $seq... *\n";
	}elsif( /x86_64/ ){
		my @tmp = split(/\s+/, $_);
		open WW, "$tmp[-1].clstr" || die "$!\n";
		while( <WW> ){
			chomp;
			if( /^>/ ){
				$num++;
				print OUT ">Cluster $num\n";
			}else{
				print OUT "$_\n";
			}
		}
		close WW;
	}
}
close OUT;
close IN;

