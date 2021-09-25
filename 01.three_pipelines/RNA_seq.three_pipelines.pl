#!/usr/bin/perl -w
use strict;

=head1 NAME

RNA_seq.three_pipelines.pl  --  the pipeline of denovo assembling RNA-seq, predicting ORFs for transcripts and qc.

=head1 DESCRIPTION

This script invokes trinity, soapdenovo_trans to assemble unigenes, cd-hit to cluster, wgd to get Ks and BUSCO to qc. 

=head1 Version

  Author: Jia Li, jiali@psb.ugent.be
  Version: 2.0  Date: 2021-4-16

=head1 Usage

    perl RNA_seq.pipe.pl [options] RNAseq.reads.list
    --trinity_para <str> set the parameter for trinity, default "--seqType fq --no_salmon --no_bowtie --max_memory 30G --CPU 6"
    --soap_para <str>    set the parameter for soapdenovo-trans, default "-K 25 -F -p 6"
    --gapcloser <str>    set the parameter for gapcloser, default "-t 6"
    --length <num>       set the minimum length for assembling, default 150    
    --cdhit_para <str>   set the parameter for cd-hit, default "-c 0.99 -T 4 -M 2000"

	--only_mixed         only run mixed sampling 

    --pipeline1          invoke the pipeline1 Tradintional Trinity method, Trinity -> Transdecoder -> cd-hit -> busco -> wgd  
    --step_p1 <num>      set the step for running pipeline1, default 12345
    --pipeline2          invoke the pipeline2 Trinity longest orf method based on the Trinity cluster, longest_orf_selectinog -> busco -> wgd
    --step_p2 <num>      set the step for running pipeline2, default 123
    --pipeline3          invoke the pipeline3 Soapdenovo_trans method, Soapdenovo_trans -> Transdecoder -> cd-hit -> busco -> wgd
    --step_p3 <num>      set the step for running pipeline3, default 12345

    --reference_genmoe   set the reference genome for rnaQUAST
    --reference_gtf      set the reference gtf for rnaQUAST
    --reference_pep      set the reference proteome for transrate

    --outdir <str>       set the result directory, default = "./"
    
    --verbose            output verbose information to screen
    --help               output help information to screen

=head1 Example

  perl RNA_seq.three_pipelines.pl RNAseq.reads.list --pipeline1 --step1 12345 --pipeline2 --step2 123 --pipeline3 --step3 12345 

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#get options from command line into variables and set default values
my ($Pipeline1, $Pipeline2, $Pipeline3);
my ($Outdir, $Step1, $Step2, $Step3);
my ($Trinity_para, $Length);
my ($Soap_para, $Gapcloser_para);
my $Cdhit_para;
my $Mixed_sampling;
my ($Ref_genome, $Ref_gtf, $Ref_pep);
my ($Verbose, $Help);

GetOptions(
	"trinity_para:s"=>\$Trinity_para,
	"length:i"=>\$Length,
	"soap_para:s"=>\$Soap_para,
	"gapcloser_para:s"=>\$Gapcloser_para,
	"cdhit_para:s"=>\$Cdhit_para,
	"outdir:s"=>\$Outdir,
	"step_p1:s"=>\$Step1,
	"step_p2:s"=>\$Step2,
	"step_p3:s"=>\$Step3,
	"pipeline1!"=>\$Pipeline1,
	"pipeline2!"=>\$Pipeline2,
	"pipeline3!"=>\$Pipeline3,
	"only_mixed!"=>\$Mixed_sampling,
	"reference_genome:s"=>\$Ref_genome,
	"reference_gtf:s"=>\$Ref_gtf,
	"reference_pep:s"=>\$Ref_pep,
	"verbose!"=>\$Verbose,
	"help!"=>\$Help
);

$Trinity_para ||= "--seqType fq --no_salmon --no_bowtie --max_memory 30G --CPU 6";
$Length ||= 150;

$Soap_para ||= "-K 25 -F -p 6";
$Gapcloser_para ||= "-t 6";


$Cdhit_para ||= "-c 0.99 -T 4 -M 2000";
my $Cdhit_C = $1 if $Cdhit_para =~ /-c\s(0\.\d+)\s.*/;
$Outdir ||= ".";
$Step1 ||= "1234567";
$Step2 ||= "34567";
$Step3 ||= "1234567";

die `pod2text $0` if (@ARGV == 0 || $Help);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my %config;
parse_config("$Bin/config.txt", \%config);
my $SOAP = $config{"soapdenovo_trans"};
my $GAPCLOSER = $config{"gapcloser"};
my $TRINITY = $config{"trinity"};
my $CDHITEST = $config{"cd_hit_est"};
my $CDHIT = $config{"cd_hit"};
my $SEQCLEAN = $config{"seqclean"};
my $TRANSRATE = $config{"transrate"};
my $RNAQUAST = $config{"rnaQUAST"};
my $SELECT_FASTA_BASED_ON_LENGTH = $config{"Select_fasta_based_on_length"};
my $BUSCO = $config{"busco"};
my $BUSCO_DB = $config{"busco_db"};
my $LONGORFS = $config{"longorfs"};
my $PREDICT = $config{"predict"};
my $DIAMOND = $config{"diamond"};
my $HMMSCAN = $config{"hmmscan"};
my $SELECTING_GFF3 = $config{"selecting_gff3"};
my $Selecting_Trinity_transcript_based_on_Transdecoder = $config{"Selecting_Trinity_transcript_based_on_Transdecoder"};
my $Selecting_Trinity_transcript_based_on_Transdecoder_longest_orfs = $config{"Selecting_Trinity_transcript_based_on_Transdecoder_longest_orfs"};
my $Select_transdecoder_cds = $config{"Select_transdecoder_cds"};
my $Select_transcripts = $config{"Select_transcripts"};
my $GET_clstr_info = $config{"get_clstr_info"};

my $pwd = `pwd`;
chomp $pwd;

my $read_lst = shift;

my (%reads, @all_reads);
open IN, $read_lst || die "No such a file: $read_lst\n";
while( <IN> ){
	chomp;
	my $fq1 = $_;
	my $fq2 = <IN>;
	chomp $fq2;
	my @tmp = split(/\//, $fq1);
	my $sp = (split(/\./, $tmp[-1]))[0];
	$reads{$sp} = [ $fq1, $fq2 ];
	push(@{$reads{"Mixed"}}, $fq1);
	push(@{$reads{"Mixed"}}, $fq2);
}
close IN;

if( $Mixed_sampling ){
	for my $i ( sort keys %reads ){
		delete $reads{$i} unless $i eq "Mixed";
	}
}


if( $Pipeline1 ){
	my %trintiy_out;
	my $wd = "$pwd/pipeline1";
	mkdir $wd unless -d $wd;
	chdir $wd;

	##### run trinity #####
	my $sub_dir1 = "$wd/01.assembling";
	my $shell1 = "$sub_dir1/S1.assembling.shell";
	run_trinity(\%reads, \%trintiy_out, $sub_dir1, $shell1);
	`sh $shell1` if ( $Step1=~ /1/);

	##### run transdecoder #####
	my $sub_dir2 = "$wd/02.ORF_finding";
	my $shell2 = "$sub_dir2/S2.transdecoder.shell";
	my %Tri_orf_output;
	run_transdecoder(\%trintiy_out, \%Tri_orf_output, $sub_dir2, $shell2);
	`sh $shell2` if $Step1=~ /2/;

	##### run cdhit #####
	my $sub_dir3 = "$wd/03.cd-hit";
	my $shell3 = "$sub_dir3/S3.cdhit.shell";
	my (%Orfs, %Cds, %Tri_select_out);
	run_cdhit(\%Tri_orf_output, \%Orfs, \%Cds, $sub_dir3, $shell3);
	`sh $shell3` if $Step1 =~ /3/;
	select_cds(\%Tri_orf_output, \%Orfs, \%trintiy_out, \%Tri_select_out);

	##### run wgd #####
	my $sub_dir4 = "$wd/04.wgd";
	my $shell4 = "$sub_dir4/S4.run_wgd.shell";
	run_wgd(\%Cds, $sub_dir4, $shell4);
	`sh $shell4` if $Step1 =~ /4/;

	##### run busco #####
	my $sub_dir5 = "$wd/05.busco";
	my $shell5 = "$sub_dir5/S5.run_busco.shell";	
	run_busco(\%Tri_select_out, \%Orfs, $sub_dir5, $shell5);
	`sh $shell5` if $Step1 =~ /5/;
	chdir $pwd;

	##### run rnaQUAST #####
	die "Please provide --reference_genome and --reference_gtf for rnaQUAST\n" unless $Ref_genome && $Ref_gtf;
	my $sub_dir6 = "$wd/06.rnaQUAST";
	my $shell6 = "$sub_dir6/S6.run_rnaquast.shell";
	my %missinfo;
	run_rnaquast(\%Tri_select_out, $Ref_genome, $Ref_gtf, \%missinfo, $sub_dir6, $shell6);
	`sh $shell6` if $Step1 =~ /6/;

	##### run transrate #####
	die "Please provide --reference_pep for transrate\n" unless $Ref_pep;
	my $sub_dir7 = "$wd/07.transrate";
	my $shell7 = "$sub_dir7/S7.run_transrate.shell";
	run_transrate(\%reads, \%Tri_select_out, $Ref_pep, $sub_dir7, $shell7);
	`sh $shell7` if $Step1 =~ /7/;
	
	##### write outfiles path #####
	open O, ">$wd/pipeline1.all_outputs.lst" || die "$!\n";
	for my $sp ( sort keys %Orfs ){
		print O "pipeline1\t$sp\t$Tri_select_out{$sp}\t$Cds{$sp}\t$Orfs{$sp}\t$wd/04.wgd/$sp/01.wgd_dmd/$sp.selected.transdecoder.cds.mcl\t$wd/04.wgd/$sp/02.wgd_ksd/$sp.selected.transdecoder.cds.ks.tsv\t$missinfo{$sp}\n";
	}
	close O;
}	


if( $Pipeline2 ){
	### Note this pipeline will use the orf_finding of pipeline1 ###
	my $pipe1_orf = "$pwd/pipeline1/02.ORF_finding";
	my $pipe1_tri_out = "$pwd/pipeline1/01.assembling";
	die "Note this pipeline will use the orf_finding output of pipeline1, please run pipeline1 first\n" unless -d $pipe1_orf && -d $pipe1_tri_out;
	my $wd = "$pwd/pipeline2";
	mkdir $wd unless -d $wd;
	chdir $wd;
	
	##### select the longest orf of each Trinity cluster #####
	my $sub_dir3 = "$wd/03.select_longest_orf_each_Trinity_cluster";
	my $shell = "$sub_dir3/S3.select.shell";
	my (%Orfs, %Cds, %Gff, %Tri_select_out);
	run_pipe2(\%reads, \%Orfs, \%Cds, \%Gff, \%Tri_select_out, $sub_dir3, $shell, $pipe1_orf);
	`sh $shell` if $Step2 =~ /3/;

	##### run wgd #####
	my $sub_dir4 = "$wd/04.wgd";
	my $shell4 = "$sub_dir4/S4.run_wgd.shell";
	run_wgd(\%Cds, $sub_dir4, $shell4);
	`sh $shell4` if $Step2 =~ /4/;

	##### run busco #####
	my $sub_dir5 = "$wd/05.busco";
	my $shell5 = "$sub_dir5/S5.run_busco.shell";	
	run_busco(\%Tri_select_out, \%Orfs, $sub_dir5, $shell5);
	`sh $shell5` if $Step2 =~ /5/;
	chdir $pwd;

	##### run rnaQUAST #####
	die "Please provide --reference_genome and --reference_gtf for rnaQUAST\n" unless $Ref_genome && $Ref_gtf;
	my $sub_dir6 = "$wd/06.rnaQUAST";
	my $shell6 = "$sub_dir6/S6.run_rnaquast.shell";
	my %missinfo;
	run_rnaquast(\%Tri_select_out, $Ref_genome, $Ref_gtf, \%missinfo, $sub_dir6, $shell6);
	`sh $shell6` if $Step2 =~ /6/;

	##### run transrate #####
	die "Please provide --reference_pep for transrate\n" unless $Ref_pep;
	my $sub_dir7 = "$wd/07.transrate";
	my $shell7 = "$sub_dir7/S7.run_transrate.shell";
	run_transrate(\%reads, \%Tri_select_out, $Ref_pep, $sub_dir7, $shell7);
	`sh $shell7` if $Step2 =~ /7/;
	
	##### write outfiles path #####
	open O, ">$wd/pipeline2.all_outputs.lst" || die "$!\n";
	for my $sp ( sort keys %Orfs ){
		print O "pipeline2\t$sp\t$Tri_select_out{$sp}\t$Cds{$sp}\t$Orfs{$sp}\t$wd/04.wgd/$sp/01.wgd_dmd/$sp.selected.transdecoder.cds.mcl\t$wd/04.wgd/$sp/02.wgd_ksd/$sp.selected.transdecoder.cds.ks.tsv\t$missinfo{$sp}\n";
	}
	close O;
}
		

if( $Pipeline3 ){
	my %soap_out;
	my $wd2 = "$pwd/pipeline3";	
	mkdir $wd2 unless -d $wd2;
	chdir $wd2;

	##### run soapdenovo_trans #####
	my $sub_dir1 = "$wd2/01.assembling";
	my $shell1 = "$sub_dir1/S1.assembling.shell";
	run_soap(\%reads, \%soap_out, $sub_dir1, $shell1);
	`sh $shell1` if $Step3 =~ /1/;

	##### run transdecoder #####
	my $sub_dir2 = "$wd2/02.ORF_finding";
	my $shell2 = "$sub_dir2/S2.transdecoder.shell";
	my %Soap_orf_output;
	run_transdecoder(\%soap_out, \%Soap_orf_output, $sub_dir2, $shell2);
	`sh $shell2` if $Step3 =~ /2/;

	##### run cdhit #####
	my $sub_dir3 = "$wd2/03.cd-hit";
	my $shell3 = "$sub_dir3/S3.cdhit.shell";
	my (%Orfs, %Cds, %Soap_select_out);
	run_cdhit(\%Soap_orf_output, \%Orfs, \%Cds, $sub_dir3, $shell3);
	`sh $shell3` if $Step3 =~ /3/;
	select_cds(\%Soap_orf_output, \%Orfs, \%soap_out, \%Soap_select_out);

	##### run wgd #####
	my $sub_dir4 = "$wd2/04.wgd";
	my $shell4 = "$sub_dir4/S4.run_wgd.shell";
	run_wgd(\%Cds, $sub_dir4, $shell4);
	`sh $shell4` if $Step3 =~ /4/;

	##### run busco #####
	my $sub_dir5 = "$wd2/05.busco";
	my $shell5 = "$sub_dir5/S5.run_busco.shell";	
	run_busco(\%Soap_select_out, \%Orfs, $sub_dir5, $shell5);
	`sh $shell5` if $Step3 =~ /5/;
	chdir $pwd;

	##### run rnaQUAST #####
	die "Please provide --reference_genome and --reference_gtf for rnaQUAST\n" unless $Ref_genome && $Ref_gtf;
	my $sub_dir6 = "$wd2/06.rnaQUAST";
	my $shell6 = "$sub_dir6/S6.run_rnaquast.shell";
	my %missinfo;
	run_rnaquast(\%Soap_select_out, $Ref_genome, $Ref_gtf, \%missinfo, $sub_dir6, $shell6);
	`sh $shell6` if $Step3 =~ /6/;

	##### run transrate #####
	die "Please provide --reference_pep for transrate\n" unless $Ref_pep;
	my $sub_dir7 = "$wd2/07.transrate";
	my $shell7 = "$sub_dir7/S7.run_transrate.shell";
	run_transrate(\%reads, \%Soap_select_out, $Ref_pep, $sub_dir7, $shell7);
	`sh $shell7` if $Step3 =~ /7/;
	
	##### write outfiles path #####
	open O, ">$wd2/pipeline3.all_outputs.lst" || die "$!\n";
	for my $sp ( sort keys %Orfs ){
		print O "pipeline3\t$sp\t$Soap_select_out{$sp}\t$Cds{$sp}\t$Orfs{$sp}\t$wd2/04.wgd/$sp/01.wgd_dmd/$sp.selected.transdecoder.cds.mcl\t$wd2/04.wgd/$sp/02.wgd_ksd/$sp.selected.transdecoder.cds.ks.tsv\t$missinfo{$sp}\n";
	}
	close O;
}



####################################
########### SUB FUNCTIONS ##########
####################################

sub parse_config{ 
	my $conifg_file = shift;
 	my $config_p = shift;

	my $error_status = 0;
	open IN, $conifg_file || die "fail open: $conifg_file";
	while( <IN> ){
		next if ( /^#/ );
        if( /(\S+)\s*=\s*(\S+)/ ){
            my ($software_name, $software_address) = ($1, $2);
            $config_p->{$software_name} = $software_address;
            if (! -e $software_address){
                warn "Non-exist: $software_name $software_address\n";
                $error_status = 1;
            }
        }
    }
    close IN;
    die "\nExit due to error of software configuration\n" if($error_status);
}


sub run_trinity{
	my $reads = shift;
	my $trinity_out = shift;
	my $dir = shift;
	my $shell = shift;
	
	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$reads ){
		my ($left, $right);
		if( $sample eq "Mixed" ){
			my @tmp_reads = @{$reads{$sample}};
			my @left_reads = @tmp_reads[ grep { $_ % 2 == 0 } 0 .. $#tmp_reads ];
			my @right_reads = @tmp_reads[ grep { $_ % 2 } 0 .. $#tmp_reads ];
			$left = join(",", @left_reads);
			$right = join(",", @right_reads);
		}else{
			$left = $reads{$sample}[0];
			$right = $reads{$sample}[1];
		}
		print SHELL "$TRINITY $Trinity_para --min_contig_length $Length --left $left --right $right --output $dir/trinity_$sample\n";
		print SHELL "sed 's/TRINITY_DN/T/' $dir/trinity_$sample/Trinity.fasta >$dir/trinity_$sample/Trinity.reID.fasta\n";
		print SHELL "perl $SELECT_FASTA_BASED_ON_LENGTH $dir/trinity_$sample/Trinity.reID.fasta $Length >$dir/trinity_$sample/$sample.trinity.fa\n";
		$trinity_out->{$sample} = "$dir/trinity_$sample/$sample.trinity.fa";
	}
	close SHELL;
}


sub run_soap{
	my $reads = shift;
	my $soap_out = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %reads ){
		my $sp_dir = "$dir/$sample";
		mkdir $sp_dir unless -d $sp_dir;
		open CONF, ">$sp_dir/$sample.soap.conf" || die "$!\n";
		print CONF "max_rd_len=150\n[LIB]\nrd_len_cutoff=150\navg_ins=200\nreverse_seq=0\nasm_flags=3\nmap_len=32\n";
		if( $sample eq "Mixed" ){
			my @tmp_reads = @{$reads{$sample}};
			my @left_reads = @tmp_reads[ grep { $_ % 2 == 0 } 0..$#tmp_reads ];
			my @right_reads = @tmp_reads[ grep { $_ % 2 } 0..$#tmp_reads ];
			for my $i ( @left_reads ){
				print CONF "q1=$i\n";
			}
			for my $i ( @right_reads ){
				print CONF "q2=$i\n";
			}
		}else{
			print CONF "q1=$reads{$sample}[0]\nq2=$reads{$sample}[1]\n";
		}
		close CONF;
		
		print SHELL "$SOAP all -s $sp_dir/$sample.soap.conf $Soap_para -o $sp_dir/$sample.soap_trans\n";
		print SHELL "$GAPCLOSER -a $sp_dir/$sample.soap_trans.scafSeq -b $sp_dir/$sample.soap.conf -o $sp_dir/$sample.soap_trans.GapCloser.fa $Gapcloser_para\n";
		print SHELL "perl $SELECT_FASTA_BASED_ON_LENGTH $sp_dir/$sample.soap_trans.GapCloser.fa $Length >$sp_dir/$sample.soap_gapcloser.le$Length.fa\n";
		$soap_out->{$sample} = "$sp_dir/$sample.soap_gapcloser.le$Length.fa";
	}
	close SHELL;
}	


sub run_transdecoder{
	my $clean_fasta = shift;
	my $orf_out = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$clean_fasta ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "ln -sf $clean_fasta->{$sample} $sub_dir/$sample.trans.clean.fa\n";
		print SHELL "cd $sub_dir; module load transdecoder\n";
		print SHELL "TransDecoder.LongOrfs -t $sub_dir/$sample.trans.clean.fa\n";
		print SHELL "TransDecoder.Predict -t $sub_dir/$sample.trans.clean.fa\n";
		$orf_out->{$sample} = "$sub_dir/$sample.trans.clean.fa.transdecoder.pep";
	}
	close SHELL;
}


sub run_cdhit{
	my $orfs = shift;
	my $cdhit_out = shift;
	my $cdhit_cds = shift;
	my $dir = shift;
	my $shell = shift;
	
	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";;
	for my $sample ( sort keys %$orfs ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "$CDHIT -i $orfs->{$sample} $Cdhit_para -o $sub_dir/$sample.transdecoder.cdhit.pep\n";
		$cdhit_out->{$sample} = "$sub_dir/$sample.transdecoder.cdhit.pep";
		$cdhit_cds->{$sample} = "$sub_dir/$sample.transdecoder.cdhit.cds";
	}
	close SHELL;
}


sub select_cds{
	my $orf_out = shift;
	my $cdhit_out = shift;
	my $ass_out = shift;
	my $ass_select_out = shift;

	for my $sp ( sort keys %$orf_out ){
		my %all_orf;
		my $all_cds = $orf_out->{$sp};
		$all_cds =~ s/pep$/cds/;
		my (%All_cds, $id);
		open IN, $all_cds || die "No such a file: $all_cds\n";
		while( <IN> ){
			chomp;
			if( />/ ){
				$id = $1 if $_ =~ />(\S+).*/;
			}else{
				$All_cds{$id} .= $_;
			}
		}
		close IN;

		my $cdhit = $cdhit_out->{$sp};
		my (%select_id, @select_ass);
		open IN, $cdhit || die "No such a file: $cdhit\n";
		while( <IN> ){
			chomp;
			if( />/ ){
				$select_id{$1} = 1 if $_ =~ />(\S+).*/;
				my $tmp_id = (split(/\./))[0];
				$tmp_id =~ s/>//;
				push(@select_ass, $tmp_id);
			}
		}
		close IN;

		my $select_out = $cdhit;
		$select_out =~ s/pep$/cds/;
		open OUT, ">$select_out";
		for my $key ( sort keys %select_id ){
			print OUT ">$key\n$All_cds{$key}\n";
		}
		close OUT;

		my %all_fasta;
		open IN, $ass_out->{$sp} || die "No such a file: $ass_out->{$sp}\n";
		while( <IN> ){
			chomp;
			if( />/ ){
				$id = $1 if $_ =~ />(\S+).*/;
			}else{
				$all_fasta{$id} .= $_;
			}
		}			
		close IN;
		
		my $select_out1 = $cdhit;
		$select_out1 =~ s/pep$/assembling_unigenes.fa/;
		open OUT, ">$select_out1";
		for my $key ( @select_ass ){
			print OUT ">$key\n$all_fasta{$key}\n";
		}
		close OUT;
		$ass_select_out->{$sp} = $select_out1;
	}
}


sub run_pipe2{
	my $reads = shift;
	my $Orfs = shift;
	my $Unigenes = shift;
	my $Gff = shift;
	my $Tri_select_out = shift;
	my $dir = shift;
	my $shell = shift;
	my $pipe1_orf = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sp ( sort keys %$reads ){
		my $tmp_dir = "$dir/$sp";
		mkdir $tmp_dir unless -d $tmp_dir;
		print SHELL "perl $Selecting_Trinity_transcript_based_on_Transdecoder_longest_orfs $pipe1_orf/$sp/$sp.trans.clean.fa $pipe1_orf/$sp/$sp.trans.clean.fa.transdecoder.pep $tmp_dir/$sp.orf_info.xls $tmp_dir/$sp.cluster.xls $tmp_dir/$sp.longest_orf.unigenes.fa $tmp_dir/$sp.longest_orf.fa\n";
		print SHELL "perl $Select_transdecoder_cds $tmp_dir/$sp.longest_orf.unigenes.fa $pipe1_orf/$sp/$sp.trans.clean.fa.transdecoder.cds >$tmp_dir/$sp.longest_orf.cds.fa\n";
		print SHELL "perl $SELECTING_GFF3 $tmp_dir/$sp.longest_orf.cds.fa $pipe1_orf/$sp/$sp.trans.clean.fa.transdecoder.gff3 >$tmp_dir/$sp.longest_orf.gff3\n";
		$Orfs->{$sp} = "$tmp_dir/$sp.longest_orf.fa";
		$Unigenes->{$sp} = "$tmp_dir/$sp.longest_orf.cds.fa";
		$Gff->{$sp} = "$tmp_dir/$sp.longest_orf.gff3";
		$Tri_select_out->{$sp} = "$tmp_dir/$sp.longest_orf.unigenes.fa";
	}
	close SHELL;
}


sub run_wgd{
	my $cds = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$cds ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "cd $sub_dir\n";
		print SHELL "awk '{print \$1}' $cds->{$sample} > $sub_dir/$sample.selected.transdecoder.cds\n";
		print SHELL "module load wgd mcl diamond mafft fasttree i-adhore; export PATH=/home/jiali/miniconda3/envs/Paml/bin:\$PATH\n";
		print SHELL "wgd dmd -I 3 $sub_dir/$sample.selected.transdecoder.cds -o $sub_dir/01.wgd_dmd --nostrictcds\n";
		print SHELL "wgd ksd $sub_dir/01.wgd_dmd/$sample.selected.transdecoder.cds.mcl $sub_dir/$sample.selected.transdecoder.cds -o $sub_dir/02.wgd_ksd\n";
		print SHELL "wgd mix -ni 100 --method bgmm $sub_dir/02.wgd_ksd/$sample.selected.transdecoder.cds.ks.tsv -o $sub_dir/04.wgd_mix\n";
	}
	close SHELL;
}


sub run_busco{
	my $Tri_select_out = shift;
	my $pep = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$Tri_select_out ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "cd $sub_dir\nln -sf $Tri_select_out->{$sample} $sub_dir/$sample.unigenes.fa\nln -sf $pep->{$sample} $sub_dir/$sample.orfs.fa\n";
		print SHELL "$BUSCO -i $sample.unigenes.fa -o unigene_$sample.busco_eukaryota_odb10 -l $BUSCO_DB/v5/eukaryota_odb10 -m tran -c 2 --offline\n";
		print SHELL "$BUSCO -i $sample.unigenes.fa -o unigene_$sample.busco_viridiplantae_odb10 -l $BUSCO_DB/v5/viridiplantae_odb10 -m tran -c 2 --offline\n";
		print SHELL "$BUSCO -i $sample.unigenes.fa -o unigene_$sample.busco_embryophyta_odb10 -l $BUSCO_DB/v5/embryophyta_odb10 -m tran -c 2 --offline\n";
		print SHELL "$BUSCO -i $sample.orfs.fa -o orf_$sample.busco_eukaryota_odb10 -l $BUSCO_DB/v5/eukaryota_odb10 -m prot -c 2 --offline\n";
		print SHELL "$BUSCO -i $sample.orfs.fa -o orf_$sample.busco_viridiplantae_odb10 -l $BUSCO_DB/v5/viridiplantae_odb10 -m prot -c 2 --offline\n";
		print SHELL "$BUSCO -i $sample.orfs.fa -o orf_$sample.busco_embryophyta_odb10 -l $BUSCO_DB/v5/embryophyta_odb10 -m prot -c 2 --offline\n";
	}
	close SHELL;
}


sub run_rnaquast{
	my $ass_select_out = shift;
	my $ref_genome = shift;
	my $ref_gtf = shift;
	my $missinfo = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$ass_select_out ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "cd $sub_dir; python3 $RNAQUAST --transcripts $ass_select_out->{$sample} --reference $ref_genome --gtf $ref_gtf -t 1 -o $sub_dir\n";
		print SHELL "grep '>' $sub_dir/$sample.*_output/$sample.*.misassembled.fasta | sed 's/>//' >$sub_dir/$sample.misassembled.id\n";
		$missinfo->{$sample} = "$sub_dir/$sample.misassembled.id";
	}
	close SHELL;
}


sub run_transrate{
	my $reads = shift;
	my $ass_select_out = shift;
	my $ref_pep = shift;
	my $dir = shift;
	my $shell = shift;

	mkdir $dir unless -d $dir;
	open SHELL, ">$shell" || die "$!\n";
	for my $sample ( sort keys %$ass_select_out ){
		my $sub_dir = "$dir/$sample";
		mkdir $sub_dir unless -d $sub_dir;
		print SHELL "$TRANSRATE --assembly $ass_select_out->{$sample} --reference $ref_pep --left $reads->{$sample}[0] --right $reads->{$sample}[1] --threads 1 --output $sub_dir\n";
	}
	close SHELL;
}
