#!/usr/bin/perl
use POSIX;
$WCHAR = 180;

return(1);

################################################################################################################################################
# This routine parses the command line (@ARGV) according to the hash array specified in @_
# variable_name => {decription=>"blah", default=>"blah", ifunreadable=>"blah", ifabsent=>"blah", store=>"blah"}
# decription will be printed if @ARGV==0
# default is the default value of $variable_name
# ifunreadable will be printed if the file $variable_name is unreadable
# ifabsent will be printed if the $variable_name is empty
# store is TRUE/FALSE value of $variable_name
# array is a flag specifying that $variable_name is an variable (undef) array (array) or hash(hash)

sub parse_command_line {
    my %hash = @_;
    foreach $key(sort keys(%hash)) {
	my %param = %{$hash{$key}};
	($trash, $param{'default'}) = split /[=\n]/, `grep ^$key= makefile` if($param{'variable'});
	$$key = $param{'default'} if($param{'default'} ne undef);
	my $obligatory = $param{'ifunreadable'} || $param{'ifabsent'};
	if(@ARGV==0 && !$param{'variable'}) {
	    print STDERR "\t-$key", ($param{'store'} ? undef : " ..."), ", ", $param{'description'};
	    print STDERR ", default=$param{'default'}" if($param{'default'} ne undef);
	    print STDERR ", obligatory" if($obligatory);
	    print STDERR ", array=$param{'array'}" if($param{'array'});
	    print STDERR "\n";
	    next;
	}
	for(my $i=0;$i<@ARGV;$i++) {
	    if($ARGV[$i] eq "-$key") {
		if($param{'store'}) {
		    $$key = $param{'store'};
		}
		else {
		    if($param{'array'}) {
		    	push @{$key}, ($param{'array'} eq "hash" ? ($ARGV[++$i], $ARGV[++$i]) : $ARGV[++$i]);
		    }
		    else {
	    	        $$key = $ARGV[++$i];
		    }
		}
	    }
	}
	print STDERR "[WARNING: $key=$$key]\n" if($$key ne $param{'default'} && !$obligatory);
	die("ERROR: $key=$$key : $param{'ifunreadable'}\n") if($param{'ifunreadable'} && ! -r $$key);
	die("ERROR: $key=$$key : $param{'ifabsent'}\n") if($param{'ifabsent'} && $$key eq undef);
    }
    exit(1) if(@ARGV==0);
}


################################################################################################################################################
# get attribute field gtf style; input = string, output = hash
sub get_attributes {
    my %res = ();
    @_[0]="@_[0];" unless(@_[0]=~/;$/);
    while(@_[0]=~/([\w\_\d]+)\s*([\"\w\:\_\,\d\-\.]*)\;/g) {
	$res{$1} = $2;
	$res{$1} =~ s/^\"(.*)\"$/$1/;
    }
    return(%res);
}

# set attribute field gtf style input = hash, output = string
sub set_attributes {
    my %hash = @_;
    my @out=();
    foreach $key(sort keys(%hash)) {
        push @out, "$key \"$hash{$key}\";";
    }
    return(join(" ", @out));
}

# same but indexfile-like, not gtf-like
sub get_features {
    my %res = ();
    while(@_[0]=~/([\w\_\d]+)\=\"{0,1}(.*?)\"{0,1}\;/g) {
	$res{$1} = $2;
	#$res{$1} =~ s/\W//g;
    }
    return(%res);
}

sub set_features {
    my %hash = @_;
    my @out=();
    foreach $key(sort keys(%hash)) {
        push @out, "$key=\"$hash{$key}\";";
    }
    return(join(" ", @out));
}

################################################################################################################################################
# This routine generates and prints a command for the makefile as follows
#
# output : input depend (script if script_required) 
# 	script before input between output after
# endpoint :: output

sub make {
    my %param=@_;
    push @{$param{'depend'}}, " $param{'script'}" if($param{'script_required'});
    $param{'script'} = "perl Perl/$param{'script'}" if($param{'script'}=~/\.pl$/);
    $param{'script'} = "Rscript R/$param{'script'}" if($param{'script'}=~/\.r$/);
    print join(" ", values(%{$param{'output'}}))," : ",join(" ", values(%{$param{'input'}}), @{$param{'depend'}}), "\n";
    print "\ttouch ", join(" ", values(%{$param{'output'}})),"\n" if($param{'touch'});
    if($param{'mkdir'}) {
	my %dirs=();
	foreach my $filename(values(%{$param{'output'}})) {
	    my @array = split /\//, $filename; pop(@array);
	    next unless(@array>0);
	    $dirs{join("/", @array)}++;
	}
	print "\tmkdir -p ", join(" ", keys(%dirs)), "\n" if(keys(%dirs)>0);
    }
    print "\t$param{'script'} ",join(" ", $param{'before'}, %{$param{'input'}}, $param{'between'}, %{$param{'output'}}, $param{'after'})," \n";
    print "$param{'endpoint'} :: ", join(" ", values(%{$param{'output'}})), "\n" if($param{'endpoint'});
    print "rm-$param{'endpoint'} ::\n\t rm -f ", join(" ", values(%{$param{'output'}})), "\n" if($param{'endpoint'});
}

sub make2 {
    my %param =@_;
    my %inputs = %{$param{'inputs'}};
    $param{'script'} = "perl Perl/$param{'script'}" if($param{'script'}=~/\.pl$/);
    foreach $key(keys(%{$param{'outputs'}})) { 
	print join(" ", values(%{$param{'outputs'}{$key}})),' '; 
    }
    print ":";
    foreach $key(keys(%{$param{'inputs'}})) { 
	print join(" ", keys(%{$param{'inputs'}{$key}})),' '; 
    }
    print "\n\t$param{'script'} $param{'before'}";
    foreach $key(keys(%{$param{'inputs'}})) {
	foreach $input(keys(%{$param{'inputs'}{$key}})) {
	    print "$key $input $param{'inputs'}{$key}{$input} ";
	}
    }
    print $param{'between'}, ' ';
    foreach $key(keys(%{$param{'outputs'}})) {
	foreach $output(keys(%{$param{'outputs'}{$key}})) {
	    print "$key $output $param{'outputs'}{$key}{$output} ";
	}
    }
    print $param{'after'},"\n";
    return unless($param{'endpoint'});
    print "$param{'endpoint'} :: ";
    foreach $key(keys(%{$param{'outputs'}})) {
	print join(" ", values(%{$param{'outputs'}{$key}})), "\n";
    }
}

################################################################################################################################################
# returns formatted fraction (inclusion)/(inclusion + exclusion) if the denomintor > minc; NA otherwise
sub frac {
    my ($INC, $EXC) = @_;
    my $TOT = $INC + $EXC;
    return($TOT > $minc ?  sprintf("%.5f", $INC/$TOT) + 0 : "NA");
}
#
################################################################################################################################################
#
# strand char (+,.,-) to integer (1, 0 -1)
sub strand_c2i {
    return 1  if(@_[0] eq "+" || @_[0] eq "1");
    return -1 if(@_[0] eq "-" || @_[0] eq "-1");
    return 0;
}
#
#strand integer (1, 0 -1) to char (+,.,-)
sub strand_i2c {
    return("+") if(@_[0] eq 1);
    return("-") if(@_[0] eq -1);
    return '.';
}
#
################################################################################################################################################

sub progressbar {
    my ($current, $last, $message) = @_;
    my $width = 2**int(log($WCHAR)/log(2));
    my $i;
    if(int(($width*($current-1))/$last) < int(($width*$current)/$last)) {
        my $k = int(($width*$current)/$last);
        print STDERR "\r$message\[";
        for($i=0;$i<$k;$i++) {print STDERR "=";}
        print STDERR ">" if($k<$width);
        for($i++;$i<$width;$i++) { print STDERR " ";}
        print STDERR "\] ",($current/$last < 0.1 ? " " : ""), int(100*$current/$last),"%";
    }
    print STDERR "\n" if($current==$last);
}

################################################################################################################################################
sub sum {
    my $s=0;
    foreach my $x(@_) {$s+=$x;}
    return($s);
}

sub avg {
    return("NA") unless(@_>0);
    return(sprintf("%.2lf", sum(@_)/@_));
}

