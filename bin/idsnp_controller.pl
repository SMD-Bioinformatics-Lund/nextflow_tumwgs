#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use JSON;

our %opt;
GetOptions( \%opt, 'vcf_sample=s', 'vcf_control=s', 'sample=s', 'control=s', 'rs_bed=s', 'help' );
my $SCORE   = -20;

my $MARKERS = 0; # originally 44 but has been cut down to 40 as off 2024-02-02 again but cound from the vcf seems more logical

print_usage(
    "ERROR: --vcf_sample, --vcf_control, --sample and --control required!\n")
  unless $opt{vcf_sample}
  and $opt{vcf_control}
  and $opt{sample}
  and $opt{control}
  and $opt{rs_bed};

unless ( -e $opt{vcf_sample} ) {
    print "ERROR: SAMPLE VCF file not found\n";
    exit;
}

if ( defined( $opt{rs_bed} ) ) {
  open( my $rsbed, '<', "$opt{rs_bed}" ) or die "Could not open file";
  while (<$rsbed>) {
    $MARKERS++;
    }
  close($rsbed);
}


my %data;

# sample vcf
if ( defined( $opt{vcf_sample} ) ) {

    if ( $opt{vcf_sample} =~ /gz$/ ) {
        open( sVCF, "zcat $opt{vcf_sample} |" )
          or die "gunzip $opt{vcf_sample}: $!";
    }
    else {
        open( sVCF, "$opt{vcf_sample}" ) or die "gunzip $opt{vcf_sample}: $!";
    }
    my %labels;
    while ( my $row = <sVCF> ) {
        chomp($row);
        my @cols = split( /\t/, $row );
        next if ( $row =~ /^##/ );
        if ( $row =~ /^#CHROM/ ) {
            for my $index ( 0 .. $#cols ) { $labels{ $cols[$index] } = $index; }
        }
        else {
            my $identifier = "".$cols[ $labels{'#CHROM'} ]."_".$cols[$labels{POS}]."_".$cols[$labels{ID}]."";
            my $id = $cols[ $labels{ID} ];
            my $gt = ( split( /:/, $cols[ $labels{FORMAT} + 1 ] ) )[0];
            my @info = split( /;/, $cols[ $labels{INFO} ] );
            my ( $reff, $refr, $altf, $altr );
            for my $c (@info) {
                chomp($c);
                if ( $c =~ /^DP4=(.*)/ ) {
                    ( $reff, $refr, $altf, $altr ) = split( /,/, $1 );
                }
            }
            $data{$identifier}{ $opt{sample} }{genotype} = $gt;
            $data{$identifier}{ $opt{sample} }{quality} =
              $cols[ $labels{QUAL} ];
            $data{$identifier}{ $opt{sample} }{dp} =
              $reff + $refr + $altf + $altr;
            $data{$identifier}{ $opt{sample} }{dp_ref} = $reff + $refr;
            $data{$identifier}{ $opt{sample} }{dp_alt} = $altf + $altr;
        }
    }
    close sVCF;
}

# control vcf
if ( defined( $opt{vcf_control} ) ) {

    unless ( -e $opt{vcf_control} ) {
        die "Error: Control VCF file '$opt{vcf_control}' does not exist.\n";
    }

    if ( $opt{vcf_control} =~ /gz$/ ) {
        open( cVCF, "zcat $opt{vcf_control} |" )
          or die "gunzip $opt{vcf_control}: $!";
    }
    else {
        open( cVCF, "$opt{vcf_control}" ) or die "gunzip $opt{vcf_control}: $!";
    }
    my %labels;
    while ( my $row = <cVCF> ) {
        chomp($row);
        my @cols = split( /\t/, $row );
        next if ( $row =~ /^##/ );
        if ( $row =~ /^#CHROM/ ) {
            for my $index ( 0 .. $#cols ) { $labels{ $cols[$index] } = $index; }
        }
        else {
            my $identifier = "".$cols[ $labels{'#CHROM'} ]."_".$cols[$labels{POS}]."_".$cols[$labels{ID}]."";
            #my $identifier = "". $cols[ $labels{'#CHROM'} ] . "_". $cols[ $labels{POS} ] . "";
            my $id = $cols[ $labels{ID} ];
            my $gt = ( split( /:/, $cols[ $labels{FORMAT} + 1 ] ) )[0];
            my @info = split( /;/, $cols[ $labels{INFO} ] );
            my ( $reff, $refr, $altf, $altr );
            for my $c (@info) {
                chomp($c);
                if ( $c =~ /^DP4=(.*)/ ) {
                    ( $reff, $refr, $altf, $altr ) = split( /,/, $1 );
                }
            }
            $data{$identifier}{ $opt{control} }{genotype} = $gt;
            $data{$identifier}{ $opt{control} }{quality} =
              $cols[ $labels{QUAL} ];
            $data{$identifier}{ $opt{control} }{dp} =
              $reff + $refr + $altf + $altr;
            $data{$identifier}{ $opt{control} }{dp_ref} = $reff + $refr;
            $data{$identifier}{ $opt{control} }{dp_alt} = $altf + $altr;
        }
    }
    close cVCF;
}

# summary table
my ( $total, $match, $mismatch ) = ( 0, 0, 0 );
my $out = "s" . $opt{sample} . "_c" . $opt{control} . ".csv";
open OUT, ">$out" or die "cannot write to outfile!";
print OUT
"VARIANT\tgt-SAMPLE\tQ-SAMPLE\tDP-SAMPLE\[ref\/alt\]\tgt-CONTROL\tQ-CONTROL\tDP-CONTROL\[ref\/alt\]\tCONCORDANCE\n";
for my $id ( sort keys %data ) {
    $total++;
    print OUT "$id\t";
    if ( $data{$id}{ $opt{sample} } ) {
        print OUT ""
          . $data{$id}{ $opt{sample} }{genotype} . "\t"
          . $data{$id}{ $opt{sample} }{quality} . "\t"
          . $data{$id}{ $opt{sample} }{dp} . '['
          . $data{$id}{ $opt{sample} }{dp_ref} . '/'
          . $data{$id}{ $opt{sample} }{dp_alt} . ']';
    }
    else { print OUT "NA\tNA\tNA"; }
    if ( $data{$id}{ $opt{control} } ) {
        print OUT "\t"
          . $data{$id}{ $opt{control} }{genotype} . "\t"
          . $data{$id}{ $opt{control} }{quality} . "\t"
          . $data{$id}{ $opt{control} }{dp} . '['
          . $data{$id}{ $opt{control} }{dp_ref} . '/'
          . $data{$id}{ $opt{control} }{dp_alt} . ']';
    }
    else { print OUT "\tNA\tNA\tNA"; }
    if (
            defined( $data{$id}{ $opt{sample} }{genotype} )
        and defined( $data{$id}{ $opt{control} }{genotype} )
        and ( $data{$id}{ $opt{sample} }{genotype} eq
            $data{$id}{ $opt{control} }{genotype} )
      )
    {
        print OUT "\tMATCH\n";
        $match++;
    }
    else { print OUT "\tMISMATCH\n"; $mismatch++; }
}

my $missing_data = $MARKERS - $total;
print OUT "\n\# missing data: $missing_data calls out of $MARKERS \("
  . sprintf( "%.2f", ( $missing_data / $MARKERS ) * 100 ) . "\)\n";
print OUT "\# mismatches: $mismatch\n";
print OUT "\# found matches: $match out of available markers\($total\) \("
  . sprintf( "%.2f", ( $match / $total ) * 100 ) . "\)\n";
print OUT "\# total matches: $match out of all markers\($MARKERS\) \("
  . sprintf( "%.2f", ( $match / $MARKERS ) * 100 ) . "\)\n";

close OUT;


my $jsonOut = "s" . $opt{sample} . "_c" . $opt{control} . ".json";
my $pp = ( $match / $MARKERS ) * 100; 
my %finalresult = (
    is_paired_sample      => JSON::true,
    pct_matching_snps     => 0 + sprintf("%.2f", $pp),
    nbr_non_matching_snps => int($mismatch),  
    nbr_matching_snps     => int($match),  
    nbr_missing_snps      => int($missing_data),  
    total_nbr_id_snps     => int($total)  
);
save_json_file (\%finalresult, $jsonOut);


sub print_usage {
    print "$_[0]\n\n" if $_[0];
    print
"USAGE: idsnp_controller.pl --vcf_sample <VCF FILE> --vcf_control <VCF FILE> --sample <ID> --control <ID> --rs_bed <rsId BED FILE>\n\n";
    print
"    --vcf_sample     FILE       Path to SAMPLE VCF file \(required\)\n\n";
    print
"    --vcf_control    FILE       Path to CONTROL VCF file \(required\)\n\n";
    print
"    --sample         ID         ID to SAMPLE VCF file \(required\)\n\n";
    print
"    --control        ID         ID to CONTROL VCF file \(required\)\n\n";
    print
"    --rs_bed        FILE       Path to bed file from rsID used for comparision \(required\)\n\n";
    print 
"    --help                      Will print this message\n\n";
    exit(0);
}    # print_usage


sub save_json_file {
    ## Save the results in the json format for dummping in the database

    my ($data_ref, $filename) = @_;
    my $json = JSON->new->allow_nonref;
    my $json_text = $json->pretty->encode($data_ref);
    open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
        print $fh $json_text;
    close $fh;
    print "Data has been written to $filename\n";
}