#!/usr/bin/perl -w
use strict;
die "perl $0 <Trinity.fa> <Transdecoder.pep> <selected_transcripts.xls> <outdir>\n" unless @ARGV == 4;


my $trinity = shift;
my $transdecoder = shift;
my $out_info = shift;
my $out_dir = shift;

my ($sequence, $trans_id, $trinity_cluster);
open IN, $trinity || die "$!\n";
while( <IN> ){
	chomp;
	if( />(.*)/ ){
		my $tt = $1;
		my @tmp = split(/\_|\s+/, $tt);
		$trinity_cluster = join("_", @tmp[0..2]);
		$trans_id = (split(/\s+/, $tt))[0];
		$sequence->{$trinity_cluster}{$trans_id} = "";
	}else{
		$sequence->{$trinity_cluster}{$trans_id} .= $_;
	}
}
close IN;


my %Decoder_group;
open TT, $transdecoder || die "$!\n";
while( <TT> ){
	chomp;
	if( />(.*)/ ){
		my $tt = $1;
		my @tmp = split(/\_|\.|\s+/, $tt);
		$trinity_cluster = join("_", @tmp[0..2]);
		$trans_id = join("_", @tmp[0..3]);
		push(@{$Decoder_group{$trinity_cluster}}, $trans_id);
	}
}
close TT;

mkdir $out_dir unless -d $out_dir;

open OUT, ">$out_info" || die "$!\n";
for my $key ( keys %Decoder_group ){
	print OUT "$key\t", join(",", @{$Decoder_group{$key}}), "\n";
	for my $i ( @{$Decoder_group{$key}} ){
		open FA, ">>$out_dir/$key.trinity_selected_seq" || die "$!\n";
		print FA ">$i\n$sequence->{$key}{$i}\n";
		close FA;
	}
}
close OUT;
