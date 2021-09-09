#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <id.fasta> <target.gff3> \n" unless @ARGV == 2;

my $id = shift;
my $gff3 = shift;

open IN, $id || die "$!\n";
my %hash;
while( <IN> ){
	chomp;
	if( />/ ){
		my @tmp = split(/\s+/, $_);
		$tmp[0] =~ s/>//;
		$hash{$tmp[0]} = 1;
	}
}
close IN;

open GFF, $gff3 || die "$!\n";
my $label;
while( <GFF> ){
	chomp;
#	print "$_\n" if /^#/;
#	print "$_\n" if /^$/;
	next if /^#/;
	next if /^$/;
	my @tmp = split(/\s+/, $_);
	if( $tmp[2] eq "gene" ){
		$label = $1 if $tmp[-1] =~ /^ID=.*~~(.*);Name=ORF.*/;
		print "$_\n" if exists $hash{$label};
	}elsif( $tmp[2] eq "mRNA" ){
		$label = $1 if $tmp[-1] =~ /^ID=(.*);Parent=.*/;
		print "$_\n" if exists $hash{$label};
	}elsif( $tmp[2] eq "exon" || $tmp[2] eq "CDS" || $tmp[2] eq "five_prime_UTR" || $tmp[2] eq "three_prime_UTR" ){
		$label = $1 if $tmp[-1] =~ /.*;Parent=(.)$/;
		print "$_\n" if exists $hash{$label};
	}
}
close GFF;

