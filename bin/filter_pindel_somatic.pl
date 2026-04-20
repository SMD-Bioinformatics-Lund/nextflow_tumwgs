#!/usr/bin/perl -w
use strict;
use warnings;

# Usage and input check
die "Usage: $0 IN_VCF OUT_VCF\n" unless @ARGV == 2;
die "$ARGV[0] does not exist or is empty!\n" unless -s $ARGV[0];

my ( $in_vcf, $out_vcf ) = @ARGV;

open( my $IN,  '<', $in_vcf )  or die "Cannot open $in_vcf: $!";
open( my $OUT, '>', $out_vcf ) or die "Cannot write to $out_vcf: $!";

while ( my $line = <$IN> ) {
    if ( $line =~ /^#/ ) {
        print $OUT $line;
        next;
    }

    chomp $line;
    my @fields = split /\t/, $line;
    #print "Processing variant at $fields[0]:$fields[1] ($fields[3] > $fields[4])\n";

    # Extract genotypes (GT) from all sample columns
    my @gts = map { ( split /:/ )[0] } @fields[ 9 .. $#fields ];
    #print "GTs: ", join(", ", @gts), "\n";

    # 1 Remove all variants if all GTs are 0/0
    my $all_ref = 1;
    foreach my $gt (@gts) {
        if ( $gt ne '0/0' ) {
            $all_ref = 0;
            print "Processing variant at $fields[0]:$fields[1] ($fields[3] > $fields[4])\n";
            print "GTs: ", join(", ", @gts), "\n";
            last;
        }
    }
    next if $all_ref;    # skip line if all are 0/0

    # 2 If exactly 2 samples: keep only if genotypes differ
    if ( @gts == 2 ) {
        my (  $normal_sample,$tumor_sample ) = @fields[ 9, 10 ];

        # Get AD field index in FORMAT
        my @format_fields = split /:/, $fields[8];    # FORMAT column
        my $ad_index = -1;
        for ( my $i = 0 ; $i < @format_fields ; $i++ ) {
            if ( $format_fields[$i] eq 'AD' ) {
                $ad_index = $i;
                last;
            }
        }

        # Skip if no AD field found
        next if $ad_index == -1;

        # Extract AD values (ref,alt) for both samples
        my @tumor_values  = split /:/, $tumor_sample;
        my @normal_values = split /:/, $normal_sample;

        my @tumor_ad  = split /,/, $tumor_values[$ad_index];
        my @normal_ad = split /,/, $normal_values[$ad_index];

        # Skip if AD parsing fails
        next unless @tumor_ad == 2 and @normal_ad == 2;

        my ( $t_ref, $t_alt ) = @tumor_ad;
        my ( $n_ref, $n_alt ) = @normal_ad;

        my $t_total = $t_ref + $t_alt;
        my $n_total = $n_ref + $n_alt;

        # Prevent division by zero
        next if $t_total <= 0;
        next if $n_total <= 10;

        print(
            "tumor_AD",  "\t", $t_total, "\t",
            "normal_AD", "\t", $n_total, "\n"
        );

        my $t_vaf = $t_alt / $t_total;
        my $n_vaf = $n_total == 0 ? 0 : $n_alt / $n_total;

        # Keep only if tumor VAF is at least 3x normal VAF
        next unless $t_vaf > 3 * $n_vaf;

        # Keep only if tumor VAF > 0.01 AND normal VAF <= 0.01
        next unless $t_vaf > 0.01 and $n_vaf <= 0.01;

    }

    # Print surviving variant
    print $OUT "$line\n";
}

close $IN;
close $OUT;
