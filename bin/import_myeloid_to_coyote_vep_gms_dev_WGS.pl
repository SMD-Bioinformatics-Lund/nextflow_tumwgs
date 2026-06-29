#!/usr/bin/perl -w
use strict;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use CMD::vcf_arr qw( parse_vcf );
use Getopt::Long;
use JSON;

my $MANE = "/data/bnf/ref/hg38/MANE.GRCh38.v0.9.summary.txt.gz";

my %opt;
GetOptions( \%opt, 'vcf=s', 'id=s', 'clarity-sample-id=s', 'clarity-pool-id=s', 'bam=s', 'group=s', 'cnv=s', 'transloc=s', 'qc=s', 'cnvprofile=s', 'build=s', 'gens=s', 'gensNorm=s' );

my( $vcf, $id ) = ( $opt{vcf}, $opt{id} );
my $genome_build = ( $opt{build} or "37" );

my @groups = split /,/, $opt{group};
#my @groups = ("mody");

# Read QC data
my @QC;
if( $opt{qc} ) {
    my @qc_files = split /,/, $opt{qc};
    foreach( @qc_files ) {
	if( -s $_ ) {
	    push @QC, read_json($_)
	}
	else {
	    print STDERR "WARNING: QC-json does not exist: $_\n.";
	}
    }
}

die "$vcf not found!" unless -s $vcf;

#################
# INSERT SAMPLE #
#################

# Connect to mongodb
my $client = MongoDB->connect();


# Prepare data to insert into sample collection
my $samples = $client->ns("coyote.samples");
my %sample_data = ( 'name'=>$id, 'groups'=>\@groups, 'time_added'=>DateTime->now, 'vcf_files'=>[$vcf],, 'genome_build'=>$genome_build );

if ( scalar @QC > 0 ) {
    $sample_data{QC} = \@QC;
}

# Add clarity information if specified
if( $opt{'clarity-sample-id'} ) {
    $sample_data{'clarity-sample-id'} =  $opt{'clarity-sample-id'};
}
if( $opt{'clarity-pool-id'} ) {
    $sample_data{'clarity-pool-id'} =  $opt{'clarity-pool-id'};
}
if( $opt{'bam'} ) {
    $sample_data{'bam'} =  $opt{'bam'};
}
if( $opt{'cnv'} ) {
    $sample_data{'cnv'} = $opt{'cnv'};
}
if( $opt{'transloc'} ) {
    $sample_data{'transloc'} = $opt{'transloc'};
}
if( $opt{'cnvprofile'} ) {
    $sample_data{'cnvprofile'} = $opt{'cnvprofile'};
}
if( $opt{'gens'} ) {
    $sample_data{'gens'} = $opt{'gens'};
}
if( $opt{'gensNorm'} ) {
    $sample_data{'gensNorm'} = $opt{'gensNorm'};
}

# Insert into collection
my $result2 = $samples->insert_one(\%sample_data);
my $SAMPLE_ID = $result2->inserted_id->value;


print STDERR "ID:".$SAMPLE_ID."\n";


###################
# Insert variants #
###################
my( $meta, $data, $sample_order  ) = parse_vcf( $vcf ); #, $id );

my $data_filtered;

# Fix/add various things to data structure
foreach( 0..scalar(@$data)-1 ) {

    # Add sample ID
    $data->[$_]->{SAMPLE_ID} = $SAMPLE_ID; #bless(\$SAMPLE_ID, "MongoDB::BSON::String");

    # Add type for pindel
    $data->[$_]->{INFO}->{TYPE} = $data->[$_]->{INFO}->{SVTYPE} if $data->[$_]->{INFO}->{SVTYPE};

    my @filters = split /;/, $data->[$_]->{FILTER};
    $data->[$_]->{FILTER} = \@filters;

    next if grep /^(FAIL_NVAF|FAIL_LONGDEL)$/, @filters;
    next if grep /^FAIL_PON_/, @filters;
    
    my @found_in = split /\|/, $data->[$_]->{INFO}->{variant_callers};
    $data->[$_]->{INFO}->{variant_callers} = \@found_in;
    
    

    my $first = 1;
#    for my $sid ( @$sample_order ) {
    for my $i ( 0.. scalar( @{ $data->[$_]->{GT} } )-1 ) {
	
	my $gt = $data->[$_]->{GT}->[$i];

        if( ( !defined($gt->{AF}) and !defined($gt->{VAF})) or !defined($gt->{DP}) or !defined($gt->{VD}) or !defined($gt->{GT}) ) {
#	    print STDERR Dumper($gt);
            die "Invalid VCF, should be aggregated with AF, DP, VD and GT";
        }

	if( defined $gt->{VAF} ) {
	    $gt->{AF} = $gt->{VAF};
	    delete $gt->{VAF};
	}
        if( $gt->{sample} =~ /^NORMAL_N/ ) {
            $gt->{sample} =~ s/^NORMAL_N//;
            $gt->{type} = "control";
        }
        elsif( $gt->{sample} =~ /^TUMOR_N/ ) {
            $gt->{sample} =~ s/^TUMOR_N//;
            $gt->{type} = "case";
        }
        else{
            $gt->{sample} =~ s/^N//;
            $gt->{type} = ( $first ? "case" : "control" );
        }
	$first = 0;
    }

    delete $data->[$_]->{vcf_str};
    delete $data->[$_]->{INFO}->{'technology.illumina'};

    push @{$data_filtered}, $data->[$_];
    
}    

#print Dumper($data);exit;


	 
 
my $variants = $client->ns("coyote.variants_idref");
my $var = $variants->with_codec( prefer_numeric => 1 );
my $result = $var->insert_many($data_filtered);
#print Dumper($result);

if( $opt{cnv} ) {
    open(CNV, $opt{cnv});
    my @cnvs;
    while(<CNV>) {
	chomp;

	my( $chr, $start, $end, $nprobes, $ratio, $strand, $genes, $panel) = split /\t/;

	# Get overlapping genes from panel, with panel info
	my %panel_info;
	my @panel_genes = split /,/, $panel;
	foreach(@panel_genes) {
	    my($gene, $class, $cnv_type ) = split /:/;
	    $panel_info{$gene} = {'gene'=>$gene, 'class'=>$class, 'cnv_type'=>$cnv_type};
	}
	
	# Get overlapping genes into an array
	my @genes = split /,/, $genes;

	my @gene_info;
	foreach my $gene (@genes) {
	    if($panel_info{$gene}) {
		push( @gene_info, $panel_info{$gene} );
	    }
	    else {
		push( @gene_info, {'gene'=>$gene} );
	    }
	}

	my %cnv = ('chr'=>$chr, 'start'=>$start, 'end'=>$end, 
		   'size'=>($end-$start), 'ratio'=>$ratio, 'genes'=>\@gene_info, 
		   'nprobes'=>$nprobes, 'SAMPLE_ID'=>$SAMPLE_ID);
	
	push( @cnvs, \%cnv);
    }

    my $variants = $client->ns("coyote.cnvs_wgs");
    my $var = $variants->with_codec( prefer_numeric => 1 );
    my $result = $var->insert_many(\@cnvs);
}


if( $opt{transloc} ) {
    my( $meta, $data, $sample_order  ) = parse_vcf( $opt{transloc} ); 
    my @filtered;

    my %mane = read_mane($MANE);
    
    # Fix/add various things to data structure
    foreach( 0..scalar(@$data)-1 ) {

	# Don't load DELs and DUPs
	if( $data->[$_]->{ALT} !~ /^</ ) {
	    # Add sample ID
	    $data->[$_]->{SAMPLE_ID} = $SAMPLE_ID;
	    $data->[$_]->{QUAL} = "" if $data->[$_]->{QUAL} eq ".";

	    my $keep_variant = 0;
	    my $mane_select;
	    my @all_new_ann;
	    my $add_mane = 0;
	    foreach my $ann (@{$data->[$_]->{INFO}->{ANN}}) {

		# Check if MANE select
		my $n_mane = 0;
		my @genes = split /&/, $ann->{Gene_ID};
		
		foreach my $ensg ( @genes ) {
		    my $enst = ($mane{$ensg}->{enst} or "NO_MANE_TRANSCRIPT");
		    $n_mane++ if $ann->{'HGVS.p'} =~ /$enst/;
		}

		my %new_ann;
		foreach my $key (keys %$ann) {
		    if( $key eq "Annotation" ) {
			foreach my $anno (@{$ann->{$key}}) {
			    $keep_variant = 1 if $anno eq "gene_fusion" or $anno eq "bidirectional_gene_fusion";
			}
		    }
		    my $key_nopoint = $key;
		    $key_nopoint =~ s/\.//g;
		    $new_ann{$key_nopoint} = $ann->{$key};
		}
		push(@all_new_ann, \%new_ann);
		if( $n_mane > 0 and $n_mane == @genes ) {
		    $mane_select = \%new_ann;
		    $add_mane = 1;
		}		
	    }
	    delete $data->[$_]->{INFO}->{ANN};
	    $data->[$_]->{INFO}->{ANN} = \@all_new_ann;
	    $data->[$_]->{INFO}->{MANE_ANN} = $mane_select if $add_mane;
	    delete $data->[$_]->{vcf_str};
	    push(@filtered, $data->[$_]) if $keep_variant;
	}
    }
#    print Dumper(\@filtered);
    
    my $transloc_coll = $client->ns("coyote.transloc");
    my $var = $transloc_coll->with_codec( prefer_numeric => 1 );
    my $result = $var->insert_many(\@filtered);

}


sub fix {
    my $str = shift;
    #$str =~ s/-/_/g;
    return $str;
}

sub read_json {
    my $fn = shift;

    print STDERR "Reading json $fn\n";

    open( JSON, $fn );
    my @json = <JSON>;
    my $decoded = decode_json( join("", @json ) );
    close JSON;

    return $decoded;
}

sub read_mane {
    my $fn = shift;
    open(MANE, "zcat $fn |");
    my %mane;
    while(<MANE>) {
	chomp;
	my @a = split /\t/;
	$a[1] =~ s/\.\d+$//;
	$a[5] =~ s/\.\d+$//;
	$a[7] =~ s/\.\d+$//;
	$mane{$a[1]} = {'enst'=>$a[7], 'refseq'=>$a[5]};
    }
    return %mane;
}
