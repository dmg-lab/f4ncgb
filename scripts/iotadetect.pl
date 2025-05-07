#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# save arguments following -h or --host in the scalar $host
# the '=s' means that an argument follows the option
# they can follow by a space or '=' ( --host=127.0.0.1 )
GetOptions( 'filter' => \(my $filter = 0 )
          , 'input=s' => \(my $input = "-")
          , 'matcher=s' => \(my $matcher = "{")
          , 'stepwidth=n' => \(my $stepwidth = 1)
          , 'verbose' => \(my $verbose = 0)
          , 'print' => \(my $print = 0)
          );

open my $file, $input or die "Could not open $input: $!";

my $run_length_sum = 0;
my $run_length_count = 0;
my $real_run_length_sum = 0;
my $real_run_length_count = 0;
my $max_run_length = 0;

while( my $line = $input eq "-" ? <STDIN> : <$file>)  {
    if(rindex($line, $matcher, 0) == 0) {
        my $start_printed = 0;
        $line = substr $line, length($matcher);

        if($verbose) {
            print "[iotadetect indices] $line";
        }
        
        my @numbers = ($line =~ /[0-9]+/g);

        # Add a large integer to the end, so the loop later always has an ending run.
        push(@numbers, ~0 >> 1);

        my $last = -1;
        my $run_length = 0;
        for my $n (@numbers) {
            if ($n > $last and $n <= $last + $stepwidth) {
                ++$run_length;
            } else {
                if ($run_length > 1) {
                    if ($print) {
                        if (!$start_printed) {
                            print "[iotadetect]";
                            $start_printed = 1;
                        }
                        print " $run_length";
                    }

                    $real_run_length_sum += $run_length;
                    ++$real_run_length_count;
                    $max_run_length = ($run_length,$max_run_length)[$run_length<$max_run_length];
                }

                $run_length_sum += $run_length;
                ++$run_length_count;
                
                $run_length = 1;
            }
            $last = $n;
        }
        if($start_printed) {
            print "\n";
        }
    } elsif (!$filter) {
        print $line;
    }
}

my $mean_run_length = sprintf "%.2f",$run_length_sum / $run_length_count;
my $mean_real_run_length = sprintf "%.2f", $real_run_length_sum / $real_run_length_count;
my $of_total = sprintf "%.2f", ($real_run_length_count / $run_length_count) * 100;

print "[iotadetect stats] of-total:$of_total % max-run-len: $max_run_length mean-len:$mean_run_length mean-real-len:$mean_real_run_length real-count:$real_run_length_count\n";

close $file;
