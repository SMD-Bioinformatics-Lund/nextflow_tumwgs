#!/usr/bin/perl -w

use strict;
use warnings;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use CMD::vcf_arr qw( parse_vcf );
use Getopt::Long;
use JSON;

# --- Configuration & Options ---

my $MANE = "/data/bnf/ref/hg38/MANE.GRCh38.v0.9.summary.txt.gz";

my %opt;
GetOptions( 
    \%opt, 
    'vcf=s', 'id=s', 'clarity-sample-id=s', 'clarity-pool-id=s', 'bam=s', 
    'group=s', 'cnv=s', 'transloc=s', 'qc=s', 'cnvprofile=s', 'build=s', 
    'gens=s', 'gensNorm=s' 
);

my $vcf          = $opt{vcf};
my $id           = $opt{id};
my $genome_build = $opt{build} || "37";
my @groups       = split /,/, ($opt{group} || '');

die "$vcf not found!" unless -s $vcf;

# --- Read QC Data ---

my @QC;
if ( $opt{qc} ) {
    my @qc_files = split /,/, $opt{qc};
    for my $file (@qc_files) {
        if ( -s $file ) {
            push @QC, read_json($file);
        } else {
            print STDERR "WARNING: QC-json does not exist: $file\n";
        }
    }
}

# ==========================================
# INSERT SAMPLE
# ==========================================

# Connect to MongoDB
my $client  = MongoDB->connect();
my $samples = $client->ns("coyote.samples");

# Prepare data to insert into sample collection
my %sample_data = ( 
    'name'         => $id, 
    'groups'       => \@groups, 
    'time_added'   => DateTime->now, 
    'vcf_files'    => [$vcf], 
    'genome_build' => $genome_build 
);

if ( @QC > 0 ) {
    $sample_data{QC} = \@QC;
}

# Dynamically add optional metadata fields if specified
my @optional_fields = qw(
    clarity-sample-id clarity-pool-id bam cnv transloc cnvprofile gens gensNorm
);
for my $field (@optional_fields) {
    if ( defined $opt{$field} ) {
        $sample_data{$field} = $opt{$field};
    }
}

# Insert into sample collection
my $sample_result = $samples->insert_one(\%sample_data);
my $SAMPLE_ID     = $sample_result->inserted_id->value;
print STDERR "ID: " . $SAMPLE_ID . "\n";


# ==========================================
# INSERT VARIANTS
# ==========================================

my ( $meta, $data, $sample_order ) = parse_vcf( $vcf );
my $data_filtered;

# Process and clean up variant data structures
for my $variant (@$data) {

    # Add sample ID reference
    $variant->{SAMPLE_ID} = $SAMPLE_ID;

    # Normalize structural variant types for pindel
    if ( $variant->{INFO}->{SVTYPE} ) {
        $variant->{INFO}->{TYPE} = $variant->{INFO}->{SVTYPE};
    }

    # Parse filters and apply exclusions
    my @filters = split /;/, $variant->{FILTER};
    $variant->{FILTER} = \@filters;

    next if grep /^(FAIL_NVAF|FAIL_LONGDEL)$/, @filters;
    next if grep /^FAIL_PON_/, @filters;
    
    # Parse variant callers
    my @found_in = split /\|/, $variant->{INFO}->{variant_callers};
    $variant->{INFO}->{variant_callers} = \@found_in;
    
    # Process genotype array
    my $first = 1;
    for my $gt ( @{ $variant->{GT} } ) {
    
        if ( (!defined($gt->{AF}) && !defined($gt->{VAF})) || !defined($gt->{DP}) || !defined($gt->{VD}) || !defined($gt->{GT}) ) {
            die "Invalid VCF, should be aggregated with AF, DP, VD and GT";
        }

        # Normalize VAF to AF
        if ( defined $gt->{VAF} ) {
            $gt->{AF} = $gt->{VAF};
            delete $gt->{VAF};
        }

        # Assign sample type designations
        if ( $gt->{sample} =~ /^NORMAL_N/ ) {
            $gt->{sample} =~ s/^NORMAL_N//;
            $gt->{type}   = "control";
        }
        elsif ( $gt->{sample} =~ /^TUMOR_N/ ) {
            $gt->{sample} =~ s/^TUMOR_N//;
            $gt->{type}   = "case";
        }
        else {
            $gt->{sample} =~ s/^N//;
            $gt->{type}   = $first ? "case" : "control";
        }
        $first = 0;
    }

    # Strip unwanted database clutter strings
    delete $variant->{vcf_str};
    delete $variant->{INFO}->{'technology.illumina'};

    push @$data_filtered, $variant;
}    
  
# Batch insert variants into MongoDB
my $variants_coll = $client->ns("coyote.variants_idref");
my $var_codec     = $variants_coll->with_codec( prefer_numeric => 1 );
$var_codec->insert_many($data_filtered);


# ==========================================
# INSERT CNVs
# ==========================================

if ( $opt{cnv} ) {
    open(my $cnv_fh, '<', $opt{cnv}) or die "Can't open CNV file $opt{cnv}: $!";
    my @cnvs;

    while ( <$cnv_fh> ) {
        chomp;
        my ( $chr, $start, $end, $nprobes, $ratio, $strand, $genes, $panel ) = split /\t/;

        # Get overlapping genes from panel with metadata
        my %panel_info;
        my @panel_genes = split /,/, ($panel || '');
        for (@panel_genes) {
            my ( $gene, $class, $cnv_type ) = split /:/;
            $panel_info{$gene} = { 'gene' => $gene, 'class' => $class, 'cnv_type' => $cnv_type };
        }
        
        # Build array of structured gene target info
        my @genes     = split /,/, ($genes || '');
        my @gene_info;
        for my $gene (@genes) {
            if ( $panel_info{$gene} ) {
                push @gene_info, $panel_info{$gene};
            } else {
                push @gene_info, { 'gene' => $gene };
            }
        }

        my %cnv = (
            'chr'       => $chr, 
            'start'     => $start, 
            'end'       => $end, 
            'size'      => ($end - $start), 
            'ratio'     => $ratio, 
            'genes'     => \@gene_info, 
            'nprobes'   => $nprobes, 
            'SAMPLE_ID' => $SAMPLE_ID
        );
        
        push @cnvs, \%cnv;
    }
    close $cnv_fh;

    my $cnv_coll  = $client->ns("coyote.cnvs_wgs");
    my $cnv_codec = $cnv_coll->with_codec( prefer_numeric => 1 );
    $cnv_codec->insert_many(\@cnvs);
}


# ==========================================
# INSERT TRANSLOCATIONS
# ==========================================

if ( $opt{transloc} ) {
    my ( $ts_meta, $ts_data, $ts_sample_order ) = parse_vcf( $opt{transloc} ); 
    my @filtered;
    my %mane = read_mane($MANE);
    
    for my $variant (@$ts_data) {

        # Process DUP,DEL and BND variants with gene fusion annotation exclusively
        if ( $variant->{INFO}->{SVTYPE} eq "DUP" || $variant->{INFO}->{SVTYPE} eq "BND" || $variant->{INFO}->{SVTYPE} eq "DEL" ) {
            
            $variant->{SAMPLE_ID} = $SAMPLE_ID;
            $variant->{QUAL}      = "" if $variant->{QUAL} eq ".";


            my $keep_variant = 0;
            my $mane_select;
            my @all_new_ann;
            my $add_mane     = 0;

            for my $ann ( @{ $variant->{INFO}->{ANN} } ) {

                my $n_mane = 0;
                my @genes  = split /&/, $ann->{Gene_ID};
                my $impact = $ann->{Annotation_Impact} // "UNKNOWN";
                my $filter = $variant->{FILTER}  // "UNKNOWN";

                next if $filter ne "PASS";
                next unless $impact eq "HIGH";
                next unless @genes == 2;
                
                for my $ensg ( @genes ) {
                    my $enst = $mane{$ensg}->{enst} || "NO_MANE_TRANSCRIPT";
                    $n_mane++ if $ann->{'HGVS.p'} =~ /$enst/;
                }

                # Construct clean annotation map, safely removing keys containing '.'
                my %new_ann;
                for my $key ( keys %$ann ) {
                    if ( $key eq "Annotation" ) {
                        for my $anno ( @{ $ann->{$key} } ) {
                            if ( $anno eq "gene_fusion" || $anno eq "bidirectional_gene_fusion" ) {
                                $keep_variant = 1;
                            }

                            #$keep_variant = 1 if $anno eq "gene_fusion" || $anno eq "bidirectional_gene_fusion";
                        }
                    }
                    my $key_nopoint = $key;
                    $key_nopoint =~ s/\.//g;
                    $new_ann{$key_nopoint} = $ann->{$key};
                }
                push @all_new_ann, \%new_ann;

                if ( $n_mane > 0 && $n_mane == @genes ) {
                    $mane_select = \%new_ann;
                    $add_mane    = 1;
                }       
            }

            delete $variant->{INFO}->{ANN};
            $variant->{INFO}->{ANN}      = \@all_new_ann;
            $variant->{INFO}->{MANE_ANN} = $mane_select if $add_mane;
            delete $variant->{vcf_str};

            push @filtered, $variant if $keep_variant;
        }
    }
    
    my $transloc_coll = $client->ns("coyote.transloc");
    my $trans_codec   = $transloc_coll->with_codec( prefer_numeric => 1 );
    $trans_codec->insert_many(\@filtered);
}


# ==========================================
# SUBROUTINES
# ==========================================

sub fix {
    my $str = shift;
    return $str;
}

sub read_json {
    my $fn = shift;
    print STDERR "Reading json $fn\n";

    open(my $json_fh, '<', $fn) or die "Can't open JSON file $fn: $!";
    my @json = <$json_fh>;
    close $json_fh;

    return decode_json( join("", @json) );
}

sub read_mane {
    my $fn = shift;

    open(my $mane_fh, '-|', "zcat $fn") or die "Can't pipe zcat $fn: $!";
    my %mane;

    while ( <$mane_fh> ) {
        chomp;
        my @a = split /\t/;
        $a[1] =~ s/\.\d+$//;
        $a[5] =~ s/\.\d+$//;
        $a[7] =~ s/\.\d+$//;
        $mane{$a[1]} = { 'enst' => $a[7], 'refseq' => $a[5] };
    }
    close $mane_fh;

    return %mane;
}