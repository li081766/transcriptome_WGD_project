#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <Trinity.fa> <Transdecoder.pep> <selected_transcripts.xls> <cluster.xls> <Selected_Trinity.fa> <Selected_Orfs.fa>\n" unless @ARGV == 6;

my $trinity = shift;
my $transdecoder = shift;
my $out_info = shift;
my $cluster_info = shift;
my $selected_trinity_out = shift;
my $selected_orf_out = shift;

my (%Sequence, $trans_id);
open IN, $trinity || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$trans_id = $1;
		$trans_id =~ s/TRINITY_DN/T/;
	}else{
		$Sequence{$trans_id} .= $_;
	}
}
close IN;

my ($Decoder, $pep_id, $len, $trinity_cluster, %Orf_sequence, $orf_id);
open TT, $transdecoder || die "$!\n";
while( <TT> ){
	chomp;
	if( />(.*)/ ){
		$len = 0;
		my $tt = $1;
		$orf_id = (split(/\s+/, $tt))[0];
		my @tmp = split(/\_|\.|\s+/, $tt);
		$trinity_cluster = join("_", @tmp[0..2]);
		$trans_id = join("_", @tmp[0..3]);
		$pep_id = (split(/\s+/, $tt))[0];
		$Decoder->{$trinity_cluster}{$trans_id} = "";
		#push(@{$Decoder_group{$trinity_cluster}}, $trans_id);
	}else{
		$Orf_sequence{$orf_id} .= $_;
		$len += length $_;
		$Decoder->{$trinity_cluster}{$trans_id} = [ $pep_id, $len ];
	}	
}
close TT;

open OUT, ">$out_info" || die "$!\n";
open CLUSTER, ">$cluster_info" || die "$!\n";
my $num_cluster = 0;
my ($longest_orf_id, $longest_orf_len, %Selected_trinity_id);
for my $cluster ( sort keys %$Decoder ){
	my $Cluster = $Decoder->{$cluster};
	print OUT ">$cluster\n";
	my $num = 0;
	for my $orf ( sort keys %$Cluster ){
		print OUT "$orf\t$Cluster->{$orf}[0]\t$Cluster->{$orf}[1]\n";
		if( $num < $Cluster->{$orf}[1] ){
			$longest_orf_id = $Cluster->{$orf}[0];
			$longest_orf_len = $Cluster->{$orf}[1];
			$num = $Cluster->{$orf}[1];
		}
	}
	print OUT "\nLongest\t$longest_orf_id\t$longest_orf_len\n\n";
	$Selected_trinity_id{$longest_orf_id} = $longest_orf_len;
	$cluster =~ s/T/TRINITY_DN/;
	my $trinity_id = (split(/\./, $longest_orf_id))[0];
	$trinity_id =~ s/T/TRINITY_DN/;
	print CLUSTER "$trinity_id\t$cluster\n";
}
close CLUSTER;
close OUT;

open FA1, ">$selected_trinity_out" || die "$!\n";
open FA2, ">$selected_orf_out" || die "$!\n";
for my $key ( sort keys %Selected_trinity_id ){
	my @tmp = split(/\_|\./, $key);
	my $trans_id = join("_", @tmp[0..3]);
	print FA1 ">$trans_id\n$Sequence{$trans_id}\n";
	print FA2 ">$key\n$Orf_sequence{$key}\n";
}
close FA2;
close FA1;
