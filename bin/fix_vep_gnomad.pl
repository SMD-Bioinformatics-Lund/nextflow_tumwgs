#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;
## tnscope, and perhaps other variantcallers create wierd multi-allelic insertions
## these variants get gnomadfrequencies per allele separated by "&" which coyote cannot handle
## For these rare occasions this script choses the most common allele and present a single AF


my $vcf = vcf2->new('file'=>$ARGV[0] );

print_header($ARGV[0]);
my %fields2fix = ("gnomADg" => 1, "gnomADg_AF" => 1, "gnomADg_AF_popmax" => 1, "gnomADg_popmax" => 1 );
my @csq_order = split(/\|/,$vcf->{meta}->{INFO}->{CSQ}->{Description});
$csq_order[0] = "Allele";
#$csq_order[@csq_order-1] =~ s/\"//;
#print join("\n",@csq_order);
while ( my $var = $vcf->next_var() ) {

    #print Dumper($var->{INFO}->{CSQ});

    
    my @csq_str = ();
    foreach my $trans ( @{ $var->{INFO}->{CSQ} }) {
        my @csq_item;
        foreach my $key ( @csq_order ) {
            if ($fields2fix{$key}) {
                my @value = split("\&",$trans->{$key});
                if (@value < 2) {
                    push @csq_item,$trans->{$key};
                    next;
                }
                my $value;
                if ($key eq 'gnomADg') {
                    $value = $value[0];
                }
                elsif ($key eq 'gnomADg_AF') {
                    $value = find_max(\@value);
                }
                elsif ($key eq 'gnomADg_AF_popmax') {
                    $value = find_max(\@value);
                }
                elsif ($key eq 'gnomADg_popmax') {
                    $value = $value[0];
                }
                push @csq_item,$value;
            }
            elsif ($key eq 'Consequence') {
                #print $key."\n";
                my $tmp = join("\&",@{$trans->{$key}});
                #print $tmp."\n";
                push @csq_item,$tmp;
            }
            else {
                push @csq_item,$trans->{$key};
            }
        }
        my $csq_item = join("\|",@csq_item);
        #print $csq_item."\n";
        push @csq_str,$csq_item;

        
    }
    $var->{INFO}->{"_CSQ_str"} = join(",",@csq_str);  
    vcfstr($var);
}


sub print_header {
    my $file = shift;

    system("zgrep ^## $file");
    #print "##FILTER=<ID=GERMLINE,Description=\"Germline variant, detected in normal sample\">\n";
    #print "##FILTER=<ID=GERMLINE_RISK,Description=\"Potential germline variant, from tumor sample\">\n";
    system("zgrep ^#CHROM $file");
	    
}

sub find_max {
    my $val = shift;
    my @val = @$val;

    my $max = 0;
    foreach my $v (@val) {
        if ($v > $max) {
            $max = $v;
        }
    }
    return $max;
}

sub vcfstr {
    my( $v ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	my $key2 = $info_key;
	$key2 = "_CSQ_str" if $info_key eq "CSQ";
	push @all_info, $info_key."=".$v->{INFO}->{$key2};
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

