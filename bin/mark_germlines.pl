#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use JSON;
use Data::Dumper;

my $MAX_VAF_NORMAL = 0.05;
my $MIN_VAF_NORMAL = 0.35;
my $MIN_VAF_TUMOR  = 0.45;
my $MIN_DP = 100;

# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'normal-id=s', 'assay=s' );
check_options( \%opt );

# Default to myeloid genes, for backwards compatibility
my %GENES;
my %assay_json;
if( $opt{'assay'}) {
    my $assay_json = read_json($opt{'assay'});
    %assay_json = %$assay_json;
    %GENES = ();
    foreach my $gene ( @{ $assay_json{"genes"} }) {
        $GENES{$gene} = 1;
    }
}

else {
    die "please provide a json defining what genes should be analyzed for GERMLINE inclusion";
}

my $vcf = vcf2->new('file'=>$opt{vcf} );

my $nid = ( $opt{'normal-id'} or 0 );
my $tid = ( $opt{'tumor-id'} or 0 );

print_header($opt{vcf});


while ( my $var = $vcf->next_var() ) {

    # Check if variant is in one of the selected genes
    my $in_relevant_gene = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if( $GENES{ $tx->{SYMBOL} } or $GENES{'ALL_GENES'} ) {
            $in_relevant_gene = 1;
            last;
        }
    }

    my $germline = 0;
    my $germline_risk = 0;
    my $not_germline = 0;

    # check whether it seems like a germline variant depending on VAF of variant
    for my $gt ( @{$var->{GT}}) {
        # If normal samples. Set GERMLINE filter if VAF > $MIN_VAF_NORMAL and in relevant gene
        if( $nid and $gt->{_sample_id} eq $nid ) {
            $germline = 1 if( $gt->{VAF} > $MIN_VAF_NORMAL and $gt->{DP} > $MIN_DP and $in_relevant_gene);
            $not_germline = 1 if $gt->{VAF} < $MAX_VAF_NORMAL and $gt->{DP} > $MIN_DP;
        }

        # If tumor samples. Set GERMLINE_RISK filter if VAF > $MIN_VAF_TUMOR
        if( $tid and $gt->{_sample_id} eq $tid ) {
            $germline_risk = 1 if( $gt->{VAF} > $MIN_VAF_TUMOR and $gt->{DP} > $MIN_DP );
        }
    }

    # if using exclusion criterias from assay_json reset germline if variant fail check
    if( $assay_json{"inclusion_score"} && $germline) {
        $germline = mini_rank($var,\%assay_json);
    }

    # Add GERMLINE filter and remove FAIL_NVAF for confirmed relevant germlines
    if ( $germline ) {
        my @new_filters = ("GERMLINE");
        foreach( split(';', $var->{FILTER}) ) {
            push @new_filters, $_ unless $_ eq "FAIL_NVAF";
        }
        $var->{FILTER} = join(";", @new_filters);
    }


    # Add GERMLINE_RISK filter to suspected germline variants in tumor only
    elsif ( $germline_risk and !$nid ) {
        my @filters = split ';', $var->{FILTER};
        push @filters, "GERMLINE_RISK";
        $var->{FILTER} = join(";", @filters);
    }

    vcfstr($var);
}


sub print_header {
    my $file = shift;

    system("zgrep ^## $file");
    print "##FILTER=<ID=GERMLINE,Description=\"Germline variant, detected in normal sample\">\n";
    print "##FILTER=<ID=GERMLINE_RISK,Description=\"Potential germline variant, from tumor sample\">\n";
    system("zgrep ^#CHROM $file");

}


sub read_json {
    my $fn = shift;

    open( JSONFILE, $fn );
    my @json = <JSONFILE>;
    my $decoded = decode_json( join("", @json ) );
    close JSONFILE;

    return $decoded;
}

sub mini_rank {
    my $var = shift;
    my $rank = shift;
    my %rank = %$rank;
    my $score = 0;

    ## check clinvar ##
    # if only one annotation handle differently, this is to support sole Benign annotations
    # otherwise find highest score acoring to rank->clinvar from JSON.
    # This let's the user define some value and let other be 0. This is ugly?
    my $clinvar = 0;
    # ClinVar comes either as a top-level INFO/CLNSIG or, when added through VEP --custom, as CSQ CLINVAR_CLNSIG
    my $clnsig_str = $var->{INFO}->{CLNSIG};
    if (!$clnsig_str) {
        for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
            if ($tx->{CLINVAR_CLNSIG}) {
                $clnsig_str = $tx->{CLINVAR_CLNSIG};
                last;
            }
        }
    }
    if ($clnsig_str && $rank{"clinvar"}) {
        my @clinsig = split(/[\/&|,]/,$clnsig_str);
        if (scalar(@clinsig) == 1) {
            if ($rank{"clinvar"}{$clinsig[0]}) {
                $clinvar = $rank{"clinvar"}{$clinsig[0]};
            }
        }
        else {
            foreach my $match (@clinsig) {
                if ($rank{"clinvar"}{$match}) {
                    if ($rank{"clinvar"}{$match} > $clinvar) {
                        $clinvar = $rank{"clinvar"}{$match};
                    }
                }
            }
        }
    }
    # add the clinvar score to score, Benign can disqualify the other two categories completely (is this ok?)
    if ($clinvar) {
        $score = $clinvar + $score;
    }
    ## check consequence ## ## check gnomad ##
    # save the most severe consequence defined in JSON
    # check if gnomad passes cutoff from JSON
    my $max_score = 0; my $gnomad = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if( $tx->{Consequence} && $rank{"consequence_cutoff"} && $rank{"consequence_score"}) {
            $max_score = conseqeunce($tx->{Consequence},\%rank);
        }
        if ( $tx->{gnomADg} && $rank{"gnomad_cutoff"} ) {
            my @afs =  split(",",$tx->{gnomADg_AF});
            if ($afs[0] < $assay_json{"gnomad_cutoff"}) {
                $gnomad = 1;
            }
        }
    }
    # add to score if most severe conseqeunce is within JSON cutoff
    if ( $assay_json{"consequence_cutoff"}) {
        if ($max_score >= $assay_json{"consequence_cutoff"}) {
            $score++;
        }
    }
    # add to score if gnomad is below cutoff
    if ($gnomad >= 1) {
        $score++;
    }

    # if score higher than or equal to inclusion score return true and retain GERMLINE status
    if ($score >= $assay_json{"inclusion_score"}) {
        return 1;
    }
    # else send back 0, and thus GERMLINE is no longer true
    else {
        return 0;
    }
}

sub conseqeunce {
    my $tx = shift;
    my $rank = shift;
    my %rank = %$rank;
    my $max_score = 0;
    foreach my $csq (@{$tx} ) {
        if ($rank{"consequence_score"}{$csq}) {
            if ( $rank{"consequence_score"}{$csq} > $max_score ) {
                $max_score = $rank{"consequence_score"}{$csq};
            }
        }
    }
    return $max_score;
}

sub check_options {
    my %opt = %{ $_[0] };

    help_text() unless $opt{vcf};

    die "File does not exist $opt{vcf}..." if ! -s $opt{vcf};

}


sub help_text {
    my $error = shift;

    print "\n\$ mark_germlines.pl --vcf INPUT_VCF --tumor-id [--normal-id}\n\n";
    print "   --vcf        Input vcf\n";
    print "   --normal-id  Normal sample ID\n";
    print "   --tumor-id   Tumor sample ID\n";
    print "   --assay      Select an assay. Determines which genes are included\n";
    print "   --assay_json Select an assay according to input json, filter uninteresing GERMLINE for select genes";
    print "\n";
    exit(0);
}


sub vcfstr {
    my( $v ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
    if( $info_key eq "CSQ" ) {
        # CSQ is parsed into an array by vcf2.pm, print the original string instead
        push @all_info, "CSQ=".( defined $v->{_csqstr} ? $v->{_csqstr} : ( $v->{INFO}->{_CSQ_str} // "" ) );
        next;
    }
    my $val = $v->{INFO}->{$info_key};
    push @all_info, $info_key."=".( defined $val ? $val : "" );
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}});


    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
    my @all_gt;
    for my $key ( @{$v->{FORMAT}} ) {
        push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
    }
    print "\t".join(":", @all_gt);
    }
    print "\n";
}
