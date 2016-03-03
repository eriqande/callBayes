use warnings;
use Bio::Cigar;
use Getopt::Std;

use Getopt::Std;
use vars qw/ %opt /;

sub init(){
	my $opt_string = 'hvdf:';
	getopts( "vhs:", \%opt ) or usage();
	usage() if $opt{h};
}


sub usage(){
print STDERR << "EOF";
    This program extracts variant sites  
    usage: $0 [-hvd] [-f file]

     -h        : this (help) message
     -v        : verbose output
     -d        : print debugging messages to stderr
     -f file   : file containing usersnames, one per line

    example: $0 -v -d -f file

    EOF
        exit;
    }
}


while (<>) {
	my @line = split "\t";
	my $cigar = Bio::Cigar->new("2M1D1M1I4M");
	print "Query length is ", $cigar->query_length, "\n";
	print "Reference length is ", $cigar->reference_length, "\n";

	my ($qpos, $op) = $cigar->rpos_to_qpos(3);
	$qpos = -1 if !$qpos;
	print "Alignment operation at reference position 3 is $op $qpos", "\n";
}

