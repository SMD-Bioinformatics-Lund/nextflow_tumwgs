#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Adds the germline calls (FILTER=GERMLINE, from mark_germlines.pl) to the somatic VCF of the same case.
#  * samples are matched by name, so the column order of the two files does not matter
#  * the germline CSQ is rewritten to the CSQ field layout of the somatic header (fields matched by
#    name, fields missing in the germline VEP run are left empty, germline-only fields are dropped)
#  * a variant present in both files is kept once, as the somatic record, with GERMLINE added to
#    FILTER and FAIL_NVAF removed (the same rule mark_germlines.pl applies)
#  * output is sorted by chromosome (natural order) and position

my %opt = ();
GetOptions( \%opt, 'somatic=s', 'germline=s' );
help_text() unless $opt{somatic} and $opt{germline};
die "File does not exist $opt{somatic}\n" unless -e $opt{somatic};
die "File does not exist $opt{germline}\n" unless -e $opt{germline};

my ( $s_hdr, $s_cols, $s_recs ) = read_vcf( $opt{somatic} );
my ( $g_hdr, $g_cols, $g_recs ) = read_vcf( $opt{germline} );

# ---- sample alignment ----
my @s_samples = @{$s_cols}[ 9 .. $#$s_cols ];
my %g_pos;
for my $i ( 9 .. $#$g_cols ) { $g_pos{ $g_cols->[$i] } = $i; }
for my $s (@s_samples) {
    die "Sample $s is in the somatic VCF but not in the germline VCF\n" unless exists $g_pos{$s};
}
die "Different number of samples in the somatic and germline VCF\n" if scalar(@s_samples) != scalar( keys %g_pos );
my @g_order = ( 0 .. 8, map { $g_pos{$_} } @s_samples );

# ---- header ----
my @out_hdr;
my %seen;
foreach my $l (@$s_hdr) {
    push @out_hdr, $l;
    my $k = hdr_key($l);
    $seen{$k} = 1 if defined $k;
}
foreach my $l (@$g_hdr) {
    my $k = hdr_key($l);
    next unless defined $k;
    next if $seen{$k};
    push @out_hdr, $l;
    $seen{$k} = 1;
}

# ---- CSQ layout ----
my $s_fmt = csq_format($s_hdr);
my $g_fmt = csq_format($g_hdr);
my @csq_map;
if ( $s_fmt and $g_fmt and join( "|", @$s_fmt ) ne join( "|", @$g_fmt ) ) {
    my %g_idx;
    for my $i ( 0 .. $#$g_fmt ) { push @{ $g_idx{ $g_fmt->[$i] } }, $i; }
    my %occ;
    for my $name (@$s_fmt) {
        my $n    = $occ{$name}++;
        my $idxs = $g_idx{$name};
        push @csq_map, ( $idxs ? ( $idxs->[$n] // $idxs->[-1] ) : -1 );
    }
}

# ---- merge ----
my %s_index;
for my $i ( 0 .. $#$s_recs ) {
    my $r = $s_recs->[$i];
    $s_index{ join( "\t", @{$r}[ 0, 1, 3, 4 ] ) } = $i;
}

my @all = @$s_recs;
my ( $n_dup, $n_new ) = ( 0, 0 );
foreach my $g (@$g_recs) {
    my @rec = @{$g}[@g_order];
    my $key = join( "\t", @rec[ 0, 1, 3, 4 ] );
    if ( exists $s_index{$key} ) {
        my $s = $all[ $s_index{$key} ];
        my @f = grep { $_ ne 'GERMLINE' and $_ ne 'FAIL_NVAF' and $_ ne '.' } split( /;/, $s->[6] );
        $s->[6] = join( ";", 'GERMLINE', @f );
        $n_dup++;
    }
    else {
        $rec[7] = remap_csq( $rec[7], \@csq_map ) if @csq_map;
        push @all, \@rec;
        $n_new++;
    }
}
print STDERR "combine_vcfs: $n_new germline records added, $n_dup already present in the somatic VCF\n";

# ---- sort ----
my @contigs = map { /^##contig=<ID=([^,>]+)/ ? $1 : () } @$s_hdr;
my %contig_rank;
$contig_rank{ $contigs[$_] } = $_ for 0 .. $#contigs;

my @keyed;
for my $i ( 0 .. $#all ) {
    my $c = $all[$i]->[0];
    my @ck = exists $contig_rank{$c} ? ( 0, $contig_rank{$c}, '' ) : chrom_key($c);
    push @keyed, [ \@ck, $all[$i]->[1], $i ];
}
my @sorted = sort {
         $a->[0][0] <=> $b->[0][0]
      || $a->[0][1] <=> $b->[0][1]
      || $a->[0][2] cmp $b->[0][2]
      || $a->[1] <=> $b->[1]
      || $a->[2] <=> $b->[2]
} @keyed;

print $_, "\n" for @out_hdr;
print join( "\t", @$s_cols ), "\n";
print join( "\t", @{ $all[ $_->[2] ] } ), "\n" for @sorted;


sub read_vcf {
    my $fn = shift;
    my $fh;
    if ( $fn =~ /\.gz$/ ) { open( $fh, "zcat $fn |" ) or die "Cannot read $fn\n"; }
    else                  { open( $fh, '<', $fn )     or die "Cannot read $fn\n"; }

    my ( @hdr, @cols, @recs );
    while ( my $line = <$fh> ) {
        chomp $line;
        next if $line =~ /^\s*$/;
        if    ( $line =~ /^##/ )    { push @hdr, $line; }
        elsif ( $line =~ /^#CHROM/ ) { @cols = split( /\t/, $line ); }
        else                        { push @recs, [ split( /\t/, $line, -1 ) ]; }
    }
    close $fh;
    die "No #CHROM line in $fn\n" unless @cols;
    return ( \@hdr, \@cols, \@recs );
}

sub hdr_key {
    my $l = shift;
    return "$1:$2" if $l =~ /^##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+)/;
    return undef;
}

sub csq_format {
    my $hdr = shift;
    foreach my $l (@$hdr) {
        return [ split( /\|/, $1 ) ] if $l =~ /^##INFO=<ID=CSQ,.*Format: ([^">]+)/;
    }
    return undef;
}

sub remap_csq {
    my ( $info, $map ) = @_;
    my @items = split( /;/, $info, -1 );
    foreach my $item (@items) {
        next unless $item =~ /^CSQ=(.*)$/s;
        my @tx = map {
            my @f = split( /\|/, $_, -1 );
            join( "|", map { $_ < 0 ? '' : ( $f[$_] // '' ) } @$map );
        } split( /,/, $1 );
        $item = "CSQ=" . join( ",", @tx );
    }
    return join( ";", @items );
}

sub chrom_key {
    my $c = shift;
    ( my $n = $c ) =~ s/^chr//i;
    return ( 1, $n, '' )  if $n =~ /^\d+$/;
    return ( 2, 0, '' )   if $n eq 'X';
    return ( 2, 1, '' )   if $n eq 'Y';
    return ( 2, 2, '' )   if $n =~ /^M/;
    return ( 3, 0, $n );
}

sub help_text {
    print "\n\$ combine_vcfs.pl --somatic SOMATIC_VCF --germline GERMLINE_VCF > combined.vcf\n\n";
    print "   --somatic    Somatic VCF (sample order and CSQ layout of this file are kept)\n";
    print "   --germline   Germline VCF, records flagged by mark_germlines.pl\n\n";
    exit(0);
}
