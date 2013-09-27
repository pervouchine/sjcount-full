#!/usr/bin/perl
use Switch;

$BAM_FREAD1 = 0x40;
$BAM_FREAD2 = 0x80;
$BAM_FREVERSE = 0x10;
$INFTY = 65535;
@STRAND = (1, -1);

#    print STDERR "This is a line-based utility for extracting splice junctions from SAM files\n";
#    print STDERR "$0 [-maxlen <max_intron_length>] [-minlen <min_intron_length>] [-margin <length>] [-read1 0/1] [-read2 0/1] [-binsize <bin_size>] [-lim <number_of_lines>]\n";
#    print STDERR "Input (STDIN) is a sam file, line by line\n";
#    print STDERR "Output (STDOUT) is tab delimited; columns are chr, beg, end, strand, overhang\n";

$binsize = $INFTY;

for(my $i=0;$i<@ARGV;$i++) {
    $read[0] = $ARGV[++$i] if($ARGV[$i] eq "-read1");
    $read[1] = $ARGV[++$i] if($ARGV[$i] eq "-read2");
    $margin  = $ARGV[++$i] if($ARGV[$i] eq "-margin");
    $maxlen  = $ARGV[++$i] if($ARGV[$i] eq "-maxlen");
    $minlen  = $ARGV[++$i] if($ARGV[$i] eq "-minlen");
    $binsize = $ARGV[++$i] if($ARGV[$i] eq "-binsize");
    $lim     = $ARGV[++$i] if($ARGV[$i] eq "-lim");
}


print STDERR "[Warning: set max intron length=$maxlen]\n" if($maxlen>0);
print STDERR "[Warning: set min intron length=$minlen]\n" if($maxlen>0);
print STDERR "[Warning: read margin set to $margin]\n" if($margin>0);
for($s=0; $s<2; $s++) {
    print STDERR "[Warning: will take reverse complement of read ", $s+1, "]\n" if($read[$s] % 2);
}
print STDERR "Warning: reads will be binned by exonic overhang, binsize=$binsize", "nt]\n" if($binsize<$INFTY);

while($line=<STDIN>) {
    ($query, $flag, $ref, $pos, $mapq, $cigar) = split /\t/, $line;
    $s = (($flag & $BAM_FREVERSE)>0);
    $strand = ($flag & $BAM_FREAD1) ? ($s + $read[0]) & 1 : ($s + $read[1]) & 1;
    @array = ();
    push @array, [$1, $2] while($cigar=~/(\d+)(\w)/g);

    $prev = $pos;
    for($i=0; $i<@array; $i++) {
	switch($array[$i]->[1]) {
	    case 'M' {	$pos += $array[$i]->[0];	}
	    case 'I' {			}
	    case 'D' {	$pos += $array[$i]->[0];	}
	    case 'N' {
			next if($i==0 || $i==@array-1);
			next unless($array[$i-1]->[1] eq "M" && $array[$i+1]->[1] eq "M");
			next if($array[$i-1]->[0] < $margin  || $array[$i+1]->[0] < $margin);
			next if($array[$i]->[0] < $minlen && $minlen > 0);
			next if($array[$i]->[0] > $maxlen && $maxlen > 0);
			$offset = int(($pos - $prev)/$binsize);
			print join("\t", $ref, $pos - 1, $pos + $array[$i]->[0], $STRAND[$strand], $offset),"\n";
			$pos += $array[$i]->[0];
			$prev = $pos;
		     }
	}
    }
    $n++;
    last if($n>$lim && $lim>0);
}
