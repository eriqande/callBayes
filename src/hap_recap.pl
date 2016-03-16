use warnings;
use strict;
use Bio::Cigar;
use Getopt::Std;

use Getopt::Std;
use vars qw/ %opt /;

sub init(){
	getopts( "hv:s:i:", \%opt ) or usage();
	usage() if $opt{h};
	print STDERR "Require path specification for VCF file: -v\n" and exit if not defined $opt{v};
	print STDERR "Require path specification for SAM file: -s\n" and exit if not defined $opt{s};
	print STDERR "Require Individual ID: -i\n" and exit if not defined $opt{i};
}


sub usage(){
print STDERR << "EOF";
    This program collect variant haplotype sites from SAM alignment file, and reports summary
 
    usage: $0 [-h] -v vcf_file -s sam_file -i int

     -h        : this (help) message
     -v file   : variant caller file - VCF format (!! assumed the position is sorted)
     -s file   : sequence alignment file - SAM format
     -i int    : individual ID (integer value) 

    example: $0 -v s1.vcf -s s1.sam -i 0

EOF
        exit;
}

init();
my $vcf;
my $hap;

open VCF, $opt{v};
while (<VCF>) {
	next if /^#/;
	my @line = split "\t";
	# for now, we only concern SNP site, (exclude, indel and complex events)

	# reference allele
	next if length($line[3]) > 1;
	# derived allele
	my $max_len = 0;
	my @snp = split ",", $line[4];
	for my $deriv (@snp) {
		$max_len = length($deriv) if length($deriv) > $max_len;
	} 
	next if $max_len > 1;

	#print $line[0], "\t", $line[3], "\t", $line[4], "\n";
	push @{$vcf->{$line[0]}}, $line[1];
	push @{$vcf->{"ref_".$line[0]}}, $line[3];
	push @{$vcf->{"der_".$line[0]}}, $snp[0];
}
close VCF;

open SAM, $opt{s};
while(<SAM>) {
	next if /^\@/;
	my @lines = split "\t";
	my $id = $lines[2];
	my $st_qpos = $lines[3];
	#skip if the alignment id is not found in the vcf
	next if not defined $vcf->{$id};

	my $cigar = Bio::Cigar->new($lines[5]);
	my @qseq = split "", $lines[9];
	my @qseq_qual = split "", $lines[10];
	
	next if $#qseq < 1;
	#print "Query length is ", $cigar->query_length, "\n";
	#print "Reference length is ", $cigar->reference_length, "\n";
	my $hapRead={};
	$hapRead->{"seq"} = "";
	my $ct=0;
	for my $rpos (@{$vcf->{$id}}) {
		my $rpos_adj = $rpos - $st_qpos +1;
		$ct++;

		if ($cigar->reference_length < $rpos_adj || $rpos_adj < 1) {
			$hapRead->{"seq"} .= "*";
			push @{$hapRead->{"qual"}}, "_";
			next;
		} 
		my ($qpos, $op) = $cigar->rpos_to_qpos($rpos_adj);

		if (not defined $qpos) {
			$hapRead->{"seq"} .= "_";
			push @{$hapRead->{"qual"}}, "_";
 		}
		else {
			$hapRead->{"seq"} .= $qseq[$qpos-1];
			push @{$hapRead->{"qual"}}, $qseq_qual[$qpos-1];
		}

		#print $qpos-1, "\t", $#qseq, "\t", $lines[0], "\t", $id, "\n";
		#print join "\t", $id, $rpos, $qseq[$qpos-1], $qseq_qual[$qpos-1], "\n" if $qpos != -1;	
	}

	$hap->{$id}->{$hapRead->{"seq"}}->{"ct"}++; 
	for my $i (0..$#{$vcf->{$id}}) {
		my $q = 10**(-(ord(${$hapRead->{"qual"}}[$i])-33)/10);
		${$hap->{$id}->{$hapRead->{"seq"}}->{"logC"}}[$i]+= log(1-$q) ; 
	 	${$hap->{$id}->{$hapRead->{"seq"}}->{"logW"}}[$i]+= log($q); 
	}

}

for my $id (keys %{$hap}){
	for my $h (keys %{$hap->{$id}}){
		print join "\t", $opt{i},
				$id, $h, $hap->{$id}->{$h}->{"ct"}, 
				(join ",", @{$hap->{$id}->{$h}->{"logC"}}),
				(join ",", @{$hap->{$id}->{$h}->{"logW"}}),
				 "\n";
	}
}

