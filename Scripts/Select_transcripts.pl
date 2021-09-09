#!/usr/bin/perl -w
use strict;
die "perl $0 <Trinity.fa> <Transdecoder.pep>\n" unless @ARGV == 2;


my $trinity = shift;
my $transdecoder = shift;

my (%sequence, $trans_id, $trinity_cluster);
open IN, $trinity || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$trans_id = $1;
	}else{
		$sequence{$trans_id} .= $_;
	}
}
close IN;


my %Decoder_group;
open TT, $transdecoder || die "$!\n";
while( <TT> ){
	chomp;
	if( />(\S+)/ ){
		my @tmp = split(/\./, $1);
		$Decoder_group{$tmp[0]} = 1;
	}
}
close TT;

for my $key ( keys %Decoder_group ){
	print ">$key\n$sequence{$key}\n";
}
