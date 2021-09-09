#!usr/bin/perl
use strict;

=head1 description

this program is to caculate the length of each chromosome of given species

=head1 example

perl get_chr_len.pl <.fa>

=cut
die `pod2text $0` if(@ARGV==0);
my $fa=shift;
my (%len, $name, $length, $row, $col, @list);
if( $fa =~ /\.gz$/ ){
	open IN, "gzip -dc $fa | " || die "$!\n";
}else{
	open IN, $fa or die "$!";
}

while(<IN>)
{
	chomp;
	if(/^>(\S+)/){
		$name=$1;
		$length=0;
		$row=0;
		push(@list, $name);
	}
	$/="(>\w+)";
	$_ =~ s/\s+//g;
	$length+=length $_ unless($row==0);
	$row++;
	$/="\n";
	$len{$name}=$length;
}
close IN;

foreach my $chr ( @list ){
	print "$chr\t$len{$chr}\n";
}
#foreach my $chr (sort keys %len){
#	print "$chr\t$len{$chr}\n";
#}
