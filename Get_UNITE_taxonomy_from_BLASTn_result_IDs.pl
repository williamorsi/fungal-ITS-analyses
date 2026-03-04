#!/usr/bin/perl -w
use strict;

my $input_file1 = $ARGV[0];
#this is the UNITE db taxonomy txt file
my $input_file2 = $ARGV[1];
#this is the BLASTn results file

my %hash_table;

open (INPUT1, $input_file1) or die "Unable to open input file, $input_file1";
open (INPUT2, $input_file2) or die "Unable to open input file, $input_file2";

while ( <INPUT1> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    
    if ($fields[0] =~ m/(\S+)/)
    {
        $hash_table{$1} = $fields[1];
    }
}

while ( <INPUT2> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    if ($fields[1] =~ m/(\S+)/)
        {
            my $UNITE_taxonomy = $1;
            my $OTUid=$fields[0];
            my $UNITE_accession = $fields[1];

            if ($hash_table{$UNITE_taxonomy})
            {
                         
                  print "$OTUid", "\t", "$UNITE_accession", "\t", "$hash_table{$UNITE_taxonomy}"; 
                    #this will print the genback accession, kog accession, and the KOG description
                }
			
            }

        }
        
close INPUT1;
close INPUT2;

