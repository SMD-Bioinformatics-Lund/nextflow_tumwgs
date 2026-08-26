#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use JSON;

# Evaluates DNAscope germline calls (already joint-genotyped via GVCFtyper,
# so GT reflects real zygosity -- no VAF-threshold heuristics needed here)
# against a per-assay gene panel, with optional ClinVar/consequence/gnomAD
# rank scoring reused from mark_germlines.pl's mini_rank. Only variants that
# are in-panel (and, when a scoring config is present, pass inclusion_score)
# are printed, tagged with INFO/GERMLINE_RANK.

my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'assay=s' );
check_options( \%opt );

my %GENES;
my %assay_json;
my $assay_json = read_json( $opt{assay} );
%assay_json = %$assay_json;
foreach my $gene ( @{ $assay_json{genes} } ) {
    $GENES{$gene} = 1;
}

my $vcf = vcf2->new( 'file' => $opt{vcf} );
my $tid = $opt{'tumor-id'};

print_header( $opt{vcf} );

while ( my $var = $vcf->next_var() ) {

    my ($tumor_gt) = grep { $_->{_sample_id} eq $tid } @{ $var->{GT} };
    next unless $tumor_gt;
    next if is_hom_ref_or_missing( $tumor_gt->{GT} );

    my $in_relevant_gene = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if ( $GENES{ $tx->{SYMBOL} } or $GENES{'ALL_GENES'} ) {
            $in_relevant_gene = 1;
            last;
        }
    }
    next unless $in_relevant_gene;

    my $score = 0;
    if ( $assay_json{inclusion_score} ) {
        $score = mini_rank( $var, \%assay_json );
        next unless $score >= $assay_json{inclusion_score};
    }

    add_info( $var, 'GERMLINE_RANK', $score );
    vcfstr($var);
}


sub is_hom_ref_or_missing {
    my $gt_str = shift;
    return 1 unless defined $gt_str;
    my @alleles = split /[\/|]/, $gt_str;
    return 1 if grep { $_ eq '.' } @alleles;
    return 1 unless grep { $_ ne '0' } @alleles;
    return 0;
}


sub mini_rank {
    my $var  = shift;
    my $rank = shift;
    my %rank = %$rank;
    my $score = 0;

    ## check clinvar ##
    my $clinvar = 0;
    if ( $rank{"clinvar"} ) {
        for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
            next unless $tx->{CLINVAR_CLNSIG};
            foreach my $match ( split( /[&|]/, $tx->{CLINVAR_CLNSIG} ) ) {
                if ( $rank{"clinvar"}{$match} and $rank{"clinvar"}{$match} > $clinvar ) {
                    $clinvar = $rank{"clinvar"}{$match};
                }
            }
        }
    }
    if ($clinvar) {
        $score = $clinvar + $score;
    }

    ## check consequence ## ## check gnomad ##
    my $max_score = 0;
    my $gnomad    = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if ( $tx->{Consequence} && $rank{"consequence_cutoff"} && $rank{"consequence_score"} ) {
            my $tx_score = conseqeunce( $tx->{Consequence}, \%rank );
            $max_score = $tx_score if $tx_score > $max_score;
        }
        if ( $tx->{gnomADg} && $rank{"gnomad_cutoff"} ) {
            my @afs = split( /[&,]/, ( $tx->{gnomADg_AF} // '' ) );
            if ( $afs[0] ne '' and $afs[0] < $rank{"gnomad_cutoff"} ) {
                $gnomad = 1;
            }
        }
    }
    if ( $rank{"consequence_cutoff"} ) {
        if ( $max_score >= $rank{"consequence_cutoff"} ) {
            $score++;
        }
    }
    if ( $gnomad >= 1 ) {
        $score++;
    }

    return $score;
}

sub conseqeunce {
    my $tx   = shift;
    my $rank = shift;
    my %rank = %$rank;
    my $max_score = 0;
    foreach my $csq ( @{$tx} ) {
        if ( $rank{"consequence_score"}{$csq} ) {
            if ( $rank{"consequence_score"}{$csq} > $max_score ) {
                $max_score = $rank{"consequence_score"}{$csq};
            }
        }
    }
    return $max_score;
}


sub print_header {
    my $file = shift;

    system("zgrep ^## $file");
    print "##INFO=<ID=GERMLINE_RANK,Number=1,Type=Integer,Description=\"Germline inclusion/rank score from germline_evaluate.pl\">\n";
    system("zgrep ^#CHROM $file");
}


sub read_json {
    my $fn = shift;
    die "please provide a json defining what genes should be evaluated for GERMLINE inclusion" unless $fn;

    open( JSONFILE, $fn ) or die "Could not open assay json $fn";
    my @json = <JSONFILE>;
    my $decoded = decode_json( join( "", @json ) );
    close JSONFILE;

    return $decoded;
}


sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}


sub check_options {
    my %opt = %{ $_[0] };

    help_text() unless $opt{vcf};
    help_text() unless $opt{'tumor-id'};

    die "File does not exist $opt{vcf}..." if ! -s $opt{vcf};
}


sub help_text {
    print "\n\$ germline_evaluate.pl --vcf INPUT_VCF --tumor-id ID --assay ASSAY_JSON\n\n";
    print "   --vcf        Input vcf (VEP-annotated, joint-genotyped DNAscope calls)\n";
    print "   --tumor-id   Tumor sample ID -- only this sample's genotype is evaluated\n";
    print "   --assay      Assay json: gene list + optional rank-scoring config\n";
    print "\n";
    exit(0);
}


sub vcfstr {
    my $v = shift;

    # NOTE: vcf2.pm parses CSQ (VEP) into an array of hashrefs (one per
    # transcript annotation), but also stashes the original, unparsed string
    # on the variant as _csqstr. We must use that raw string here instead of
    # reconstructing from the parsed structure, since re-deriving the
    # pipe-delimited field order from a hash isn't reliable.
    my @all_info;
    for my $info_key ( @{ $v->{INFO_order} } ) {
        if ( $info_key eq "CSQ" and defined $v->{_csqstr} ) {
            push @all_info, "CSQ=" . $v->{_csqstr};
        }
        else {
            my $val = $v->{INFO}->{$info_key};
            $val = "" if !defined($val);
            push @all_info, $info_key . "=" . $val;
        }
    }

    my @fields = (
        $v->{CHROM}, $v->{POS}, $v->{ID}, $v->{REF}, $v->{ALT}, $v->{QUAL}, $v->{FILTER},
        join( ";", @all_info ),
        join( ":", @{ $v->{FORMAT} } ),
    );

    for my $gt ( @{ $v->{GT} } ) {
        my @all_gt;
        for my $key ( @{ $v->{FORMAT} } ) {
            push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "" );
        }
        push @fields, join( ":", @all_gt );
    }

    print join( "\t", @fields ) . "\n";
}
