#!/usr/bin/perl -w
use strict;
die "\n\tUsage: perl <accession_id.lst>\n\n" unless @ARGV == 1;

my $in = shift;

########################
### accession_id.lst ###
### stem  SRR3993761 ###
### leaf  SRR3993754 ###
### root  SRR3993762 ###
########################

my $pwd = `pwd`;
chomp $pwd;
my $tmp = "$pwd/tmp";
mkdir $tmp unless -d $tmp;

my %access;
open IN, $in || die "No such a file: $in\n";
while( <IN> ){
	chomp;
	my @tmp = split(/\s+/, $_);
	$access{$tmp[0]} = $tmp[1];
}
close IN;

open SH, ">Download_Clean.shell" || die "$!\n";
for my $sp ( sort keys %access ){
	print SH "parallel-fastq-dump --sra-id $access{$sp} --tmpdir $tmp --threads 1 --outdir $pwd --split-files --gzip\n";
	print SH "trimmomatic PE -summary $pwd/$sp.$access{$sp}.trimmomatic.summary.log $pwd/$access{$sp}_1.fastq.gz $pwd/$access{$sp}_2.fastq.gz $pwd/$sp.$access{$sp}.trim_clean.1.fq.gz $pwd/$sp.$access{$sp}.trim_unpaired.1.fq.gz $pwd/$sp.$access{$sp}.trim_clean.2.fq.gz $pwd/$sp.$access{$sp}.trim_unpaired.2.fq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:50\n";
}
close SH;

`sh Download_Clean.shell`; 
