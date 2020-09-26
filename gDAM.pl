#!/usr/bin/perl

#
#	gDAM.pl     version 1.01
#
#   heroen.verbruggen@gmail.com
#
# 
#   Version history:
#   1.01  initial release

use strict;
use warnings;
use SVG;
use Graphics::ColorUtils;


### global variables #################################################################################################################################################

my (
	$infile,
	$outfile,
	$usecolor,
	$usegradient,
	$rowheight,
	
	$seqs,
	$seq_order,
	$seq_length,
	$nr_seqs,
	$charsets,
	
	$stats,
);

my $gradient_type = 'royg';
my $values = {'gradmin' => 0, 'gradmax' => 1};


### parsing command line arguments ###############################################################################################################################

unless (($ARGV[0]) and (substr($ARGV[0],0,1) eq "-") and
		($ARGV[1]) and
		($ARGV[2]) and (substr($ARGV[2],0,1) eq "-") and
		($ARGV[3]))
	{usage()}

for (my $i=0; $i<scalar(@ARGV); $i+=2) {
	if    ($ARGV[$i] eq "-i") {$infile = $ARGV[$i+1]}
	elsif ($ARGV[$i] eq "-o") {$outfile = $ARGV[$i+1]}
	elsif ($ARGV[$i] eq "-c") {$usecolor = $ARGV[$i+1]}
	elsif ($ARGV[$i] eq "-g") {$usegradient = $ARGV[$i+1]}
	elsif ($ARGV[$i] eq "-rh") {$rowheight = $ARGV[$i+1]}
	else {usage()}
}

unless ($infile and $outfile) {usage();}
unless (-e $infile) {die "\n#### FATAL ERROR ####\nfile not found: $infile\n"}
unless (defined $usecolor) {$usecolor = 1}
unless (defined $usegradient) {$usegradient = 1}
unless (defined $rowheight) {$rowheight = 150}

sub usage {
	print "\nusage:\n"; 
	print "\nmandatory parameters\n";
	print "   -i    input alignment (nexus format)\n";
	print "   -o    output file (SVG format)\n";
	print "\noptional parameters\n";
	print "   -c    1 or 0; use color (1) or gray (0)\n";
	print "   -g    1 or 0; indicate percentage site presence with color gradient\n";
	print "   -rh   row heigth (default: 150)\n";
	print "           modifying this allows you to control the height-width ratio of the table\n";
	print "\n";
	exit;
}


### parsing input alignment file ###############################################################################################################################

($seq_order,$seqs) = read_nexus_data($infile);
unless ($seq_order) {die "#### FATAL ERROR ####\nno sequences could be parsed from file\n";}
$charsets = check_charsets(read_charsets($infile));
$stats = calculate_statistics();

### generating SVG ###############################################################################################################################################

print "\ngenerating SVG ...\n";
my $xoffset = int($rowheight * get_maximal_taxon_name_length() / 3.5) + ($rowheight / 2);
my $svg = SVG->new('width' => $xoffset + $seq_length + 10 , 'height' => $rowheight * (4 + scalar(@$charsets)));
for (my $taxon_nr = 0; $taxon_nr < scalar @$seq_order; ++$taxon_nr) {
	my $taxon = $seq_order->[$taxon_nr];
	foreach my $charset (@$charsets) {
		my $seq = get_sequence($taxon,$charset);
		my ($start,$end) = determine_start_end($seq);
		my $length = $end - $start;
		unless ($seq =~ /^\-*$/) {
			if ($usecolor) {
				$svg->rect(   # draw data presence rectangle for this cell, in color
					'x' => $xoffset + $charset->{'begin'} + $start,
					'y' => $rowheight * ($taxon_nr + 1),
					'width' => $length,
					'height' => $rowheight,
					'style' => "fill:rgb(127,178,255);stroke:rgb(0,0,0);stroke-width:0.5",
				);
			} else {
				$svg->rect(  # draw data presence rectangle for this cell, in gray
					'x' => $xoffset + $charset->{'begin'} + $start,
					'y' => $rowheight * ($taxon_nr + 1),
					'width' => $length,
					'height' => $rowheight,
					'style' => "fill:rgb(100,100,100);stroke:rgb(0,0,0);stroke-width:0.5",
				);
			}
			$svg->text(  # print number of characters in this cell
				'x' => $xoffset + $charset->{'begin'} + $start + ($rowheight / 4),
				'y' => $rowheight * ($taxon_nr + 1.65),
				'-cdata' => $length,
				'font-size'	=> $rowheight/2,
				'fill' => 'rgb(255,255,255)'
			);
		}
	}
}
for (my $taxon_nr = 0; $taxon_nr < scalar @$seq_order; ++$taxon_nr) {
	my $taxon = $seq_order->[$taxon_nr];
	my $val = $stats->{'taxa'}->{$taxon}->{'sites'};
	if ($usegradient) {   # if desired, the row header gets a background color
		my $c = calc_color(1 - $val);
		my $charset = $charsets->[0];
		my $begin = $charset->{'begin'} || 0;
		$svg->rect(
			'x' => 0,
			'y' => $rowheight * ($taxon_nr + 1),
			'width' => $xoffset + $begin,
			'height' => $rowheight,
			'style' => "fill:rgb(".$c->[0].",".$c->[1].",".$c->[2].");stroke:rgb(0,0,0);stroke-width:0",
		);
	}
	$svg->text(  # print taxon name in row header
		'x' => $rowheight / 2,
		'y' => $rowheight * ($taxon_nr + 1.65),
		'-cdata' => $taxon,
		'font-size'	=> $rowheight/2,
	);
	$svg->text(  # print percentage of sites present in taxon set to the right of the row
		'x' => $xoffset + $seq_length + 1 + ($rowheight / 2),
		'y' => $rowheight * ($taxon_nr + 1.65),
		'-cdata' => sprintf("%.0f",(100 * $val)).'% of sites',
		'font-size'	=> $rowheight/2,
	);
}
for (my $taxon_nr = 0; $taxon_nr < scalar @$seq_order; ++$taxon_nr) {
	$svg->line(  # horizontal line to delimit rows
			'x1' => 0,
			'y1' => $rowheight * ($taxon_nr + 2),
			'x2' => $xoffset + $seq_length + 1 + ($rowheight / 2) + int($rowheight * 13 / 3.5),
			'y2' => $rowheight * ($taxon_nr + 2),
			'style' => 'stroke:black;stroke-width:0.5',
	);
}
$svg->line(  # horizontal line to delimit rows
		'x1' => 0,
		'y1' => $rowheight,
		'x2' => $xoffset + $seq_length + 1,
		'y2' => $rowheight,
		'style' => 'stroke:black;stroke-width:0.5',
);
foreach my $charset (@$charsets) {
	my $val = $stats->{'genes'}->{$charset}->{'sites'};
	if ($usegradient) {   # if desired, the column header gets a background color
		my $c = calc_color(1 - $val);
		$svg->rect(
			'x' => $xoffset + $charset->{'begin'},
			'y' => 0,
			'width' => $charset->{'length'},
			'height' => $rowheight,
			'style' => "fill:rgb(".$c->[0].",".$c->[1].",".$c->[2].");stroke:rgb(0,0,0);stroke-width:0",
		);
	}
	$svg->text(  # print name of character set in column header
		'x' => $xoffset + $charset->{'begin'} + ($rowheight / 2),
		'y' => $rowheight * 0.65,
		'-cdata' => $charset->{'name'}.' ('.$charset->{'length'}.')',
		'font-size'	=> $rowheight/2,
	);
	$svg->text(  # print percentage of sites present in character set underneath column
		'x' => $xoffset + $charset->{'begin'} + ($rowheight / 2),
		'y' => $rowheight * (1.65 + scalar(@$seq_order)),
		'-cdata' => sprintf("%.0f",(100 * $val)).'% of sites',
		'font-size'	=> $rowheight/2,
	);
}
foreach my $charset (@$charsets) {  # vertical line to delimit columns
	$svg->line(
		'x1' => $xoffset + $charset->{'begin'},
		'y1' => 0,
		'x2' => $xoffset + $charset->{'begin'},
		'y2' => $rowheight * (2 + scalar(@$seq_order)),
		'style' => 'stroke:black;stroke-width:0.5',
	);
}
$svg->line(  # vertical line to delimit columns
	'x1' => $xoffset + $seq_length + 1,
	'y1' => 0,
	'x2' => $xoffset + $seq_length + 1,
	'y2' => $rowheight * (2 + scalar(@$seq_order)),
	'style' => 'stroke:black;stroke-width:0.5',
);
if ($usegradient) {  # draw color gradient legend
	my $gradient = $svg->gradient(
		'-type' => "linear",
		'id' => 'leggrad',
		'x1' => "0%",
		'y1' => "0%",
		'x2' => "100%",
		'y2' => "0%",
	);
	for (my $i=0; $i<=100; $i += 10) {  # create 10 steps in gradient
		my $val = $values->{'gradmin'} + ( ($i/100) * ($values->{'gradmax'} - $values->{'gradmin'}) );
		my $c = calc_color(1 - $val);
		$gradient->stop(
			'offset' => $i."%",
			'style' => 'stop-color:rgb('.$c->[0].",".$c->[1].",".$c->[2].');stop-opacity:1',
		);
	}
	$svg->rect(  # gradient legend rectangle
		'x'			=>	($rowheight / 2),
		'y'			=>	$rowheight * (3 + scalar(@$seq_order)),
		'width'		=>	$rowheight * 5,
		'height'	=>	$rowheight / 2,
		'style' 	=> 'fill:url(#leggrad);stroke:rgb(0,0,0)',
		'stroke-width'	=>	0.5,
	);
	for (my $i=0; $i<=1; $i += 0.25) {  # draw 5 tickmarks with labels
		$svg->line(   # tickline
			'x1'	=>	($rowheight / 2) + ($rowheight * 5 * $i),
			'y1'	=>	($rowheight * (3 + scalar(@$seq_order))) - (0.1 * $rowheight),
			'x2'	=>	($rowheight / 2) + ($rowheight * 5 * $i),
			'y2'	=>	$rowheight * (3 + scalar(@$seq_order)),
			'style'	=>	'stroke:black;stroke-width:0.5',
		);
		$svg->text(   # label: percentage
			'x' => ($rowheight / 2) + ($rowheight * 5 * $i) - (int(0.3 * $rowheight * length(sprintf("%.0f",(100 * $i)).'%') / 3.5)),
			'y' => ($rowheight * (3 + scalar(@$seq_order))) - (0.15 * $rowheight),
			'-cdata' => sprintf("%.0f",(100 * $i)).'%',
			'font-size'	=> $rowheight/3,
		);
	}
	$svg->text(  # legend caption
		'x' => ($rowheight / 2),
		'y' => ($rowheight * (3 + scalar(@$seq_order))) - (0.7 * $rowheight),
		'-cdata' => 'Percentage of sites present',
		'font-size'	=> $rowheight/2,
	);
}

open FH,">$outfile" || die "\n#### FATAL ERROR ####\ncannot write to $outfile\n";
print FH $svg->xmlify();
close FH;
print "  saved to $outfile\n";

print "\n";
exit;



sub get_sequence {
	my $taxon = shift;
	my $charset = shift;
	my $seq = $seqs->{$taxon};
	return substr($seq,$charset->{'begin'}-1,$charset->{'length'})
}

sub determine_start_end {
	my $seq = shift;
	if ($seq =~ /^\-*$/) {return (0,0)}
	$seq =~ /^(\-*).*(\-*)$/;
	my ($start,$end) = ($1,$2);
	unless ($start) {$start = ''}
	unless ($end) {$end = ''}
	return (length $start, ((length $seq) - (length $end)));
}

sub get_maximal_taxon_name_length {
	my $max; $max = 0;
	foreach my $taxon (@$seq_order) {
		if (length $taxon > $max) {$max = length $taxon}
	}
	return $max;
}

sub read_nexus_data {
	my $infile = shift;
	my ($sequence_order,$seqs);
	print "\nreading alignment ...\n";
	open(FH,$infile) || die "\n#### FATAL ERROR ####\nunable to open $infile";
	my @filedump = <FH>; close FH;
	my $filestring = join('',@filedump);
	unless ($filestring =~ /begin data\;.+?matrix\s*([^\;]+)/si) {die "\n#### FATAL ERROR ####\ncould not find data block in $infile"}
	my @a = split /\n/,$1;
	foreach my $line (@a) {
		my @b = split /\s+/,$line;
		if (scalar @b > 2) {die "\n#### FATAL ERROR ####\nproblem with input: you probably have spaces in your taxon names or sequences";}
		if (scalar @b == 2) {
			unless ($seqs->{$b[0]}) {
				push @$sequence_order,$b[0];
			}
			$b[1] =~ s/\s//g;
			$seqs->{$b[0]} .= $b[1];
		}
	}
	my @keys = keys %$seqs;
	$seq_length = length($seqs->{$keys[0]});
	foreach my $id (keys %$seqs) {
		unless ($seq_length == length($seqs->{$id})) {die "\n#### FATAL ERROR ####\nsequences not of equal length"}
	}
	$nr_seqs = scalar keys %$seqs;
	print "  $nr_seqs sequences\n  $seq_length characters\n";
	return ($sequence_order,$seqs);
}

sub check_charsets {
	my $charsets = shift;
	unless ($charsets) {  # if charsets not user-defined, assume single charset
		print "\n#### WARNING ####\nno character sets were found: assuming a single character set\n";
		return [ {
			'name' => 'undefined character set',
			'begin' => '1',
			'end' => length($seqs->{$seq_order->[0]}),
			'length' => length($seqs->{$seq_order->[0]}),
		} ];
	}
	# code below dies if charsets overlap or if not all sites are accounted for
	my $sites;
	foreach my $charset (@$charsets) {
		for (my $i = $charset->{'begin'}; $i <= $charset->{'end'}; ++$i) {
			if ($sites->{'s'.$i}) {
				die "\n#### FATAL ERROR ####\nthe character sets overlap\n"
			} else {
				$sites->{'s'.$i} = 1
			}
		}
	}
	unless ($seq_length == scalar keys %$sites) {die "\n#### FATAL ERROR ####\nnot all sites are assigned to a character set\n"}
	return $charsets;
}

sub read_charsets {
	my $infile = shift;
	my $charsets;
	print "\nreading character sets ...\n";
	open(FH,$infile) || die "\n#### FATAL ERROR ####\nunable to open $infile";
	my @filedump = <FH>; close FH;
	my $filestring = join('',@filedump);
	$filestring =~ s/[\n\r]//g;
	while ($filestring =~ /charset (.+?);/ig) {
		$1 =~ /^\s*(.+?)\s*\=\s*(.+?)\s*$/;
		my ($name,$def) = ($1,$2);
		unless ($def =~ /(\d+)\-(\d+)/) {die "\n#### FATAL ERROR ####\nonly character sets of format xxx-yyy are allowed in this program\n  e.g. charset gene1 = 1-462\n"}
		my ($begin,$end) = ($1,$2);
		push @$charsets, {'name' => $name, 'begin' => $begin, 'end' => $end, 'length' => $end - $begin + 1};
		print "  charset ",scalar(@$charsets),": $name\n";
	}
	return $charsets;
}

sub calculate_statistics {
	print "\ncalculating data availability statistics\n";
	my $stats;
	$stats->{'overall'}->{'presence'} = 0;
	$stats->{'overall'}->{'sites'} = 0;
	foreach my $charset (@$charsets) {
		$stats->{'genes'}->{$charset}->{'presence'} = 0;
		$stats->{'genes'}->{$charset}->{'sites'} = 0;
	}
	for (my $taxon_nr = 0; $taxon_nr < scalar @$seq_order; ++$taxon_nr) {
		my $taxon = $seq_order->[$taxon_nr];
		$stats->{'taxa'}->{$taxon}->{'presence'} = 0;
		$stats->{'taxa'}->{$taxon}->{'sites'} = 0;
		foreach my $charset (@$charsets) {
			my $seq = get_sequence($taxon,$charset);
			unless ($seq =~ /^\-*$/) {
				my ($start,$end) = determine_start_end($seq);
				my $length = $end - $start;
				$stats->{'taxa'}->{$taxon}->{'presence'} += 1 / (scalar @$charsets);
				$stats->{'taxa'}->{$taxon}->{'sites'} += $length / $seq_length;
				$stats->{'genes'}->{$charset}->{'presence'} += 1 / $nr_seqs;
				$stats->{'genes'}->{$charset}->{'sites'} += $length / ($charset->{'length'} * $nr_seqs);
				$stats->{'overall'}->{'presence'} += 1 / ($nr_seqs * (scalar @$charsets)) ;
				$stats->{'overall'}->{'sites'} += $length;
			}
		}
	}
	$stats->{'overall'}->{'sites'} /= $nr_seqs * $seq_length;
	print "  ",sprintf("%.2f",(100 * $stats->{'overall'}->{'presence'})),"% in gene x taxon context\n";
	print "  ",sprintf("%.2f",(100 * $stats->{'overall'}->{'sites'})),"% in site x taxon context\n";
	return $stats;
}


sub calc_color {
	my $val = shift;
	my $out;
	my ($min,$max) = ($values->{'gradmin'},$values->{'gradmax'});
	if ($gradient_type =~ /^bw$/i) {
		my $r = (255 / ($min - $max));
		my $a = (255 * $max) / ($max - $min);
		my $y = int(($r * $val) + $a);
		$out = [$y,$y,$y];
	} elsif ($gradient_type =~ /^royg$/i) {
		my $c = (120 / ($min - $max));
		my $a = (120 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^rainbow$/i) {
		my $c = (300 / ($min - $max));
		my $a = (300 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^roygb$/i) {
		my $c = (200 / ($min - $max));
		my $a = (200 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^maxent$/i) {
		my $c = (240 / ($min - $max));
		my $a = (240 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^blg$/i) {
		my $a = (129 - 244) / ($max - $min);
		my $k = 244 - ($a * $min);
		my $y = ($a * $val) + $k;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^bly$/i) {
		my $r = int((255 / ($max-$min)) * ($val - $min));
		my $g = int(((155 / ($max-$min)) * ($val - $min)) + 100);
		my $b = int((255 / ($min-$max)) * ($val - $max));
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^blw$/i) {
		my $r = int((((255-20)*$val)/($min-$max)) + 20 - (((255-20)*$max)/($min-$max)));
		my $g = int((((255-0)*$val)/($min-$max)) + 0 - (((255-0)*$max)/($min-$max)));
		my $b = int((((255-215)*$val)/($min-$max)) + 215 - (((255-215)*$max)/($min-$max)));
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^rw$/i) {
		my $a = 255 / ($min - $max);
		my $c = 0 - ($a * $max);
		my $r = 255;
		my $g = int( ($a * $val) + $c );
		my $b = int( ($a * $val) + $c );
		$out = [$r,$g,$b];
	} elsif ($gradient_type =~ /^g\[(.+?)\-(.+?)\]$/) {
		my ($low,$high) = ($1,$2);
		unless (($low <= 1) and ($low >= 0) and ($high <= 1) and ($high >= 0)) {fatalerror("Values given for gradient out of 0-1 range\n");}
		my $a = 255 * (1-$low);
		my $b = 255 * (1-$high);
		my $y = int($a + ((($b - $a)/($max - $min))*($val - $min)));
		$out = [$y,$y,$y];
	} elsif ($gradient_type =~ /^ocsst$/i) {  # ocean colors sea surface temperature approximation
		my ($r,$g,$b);
		my $t = $val;
		if ($t<1) {$r=-46*$t + 54} elsif ($t<20) {$r=8} elsif ($t<27) {$r=-0.0956*$t*$t*$t*$t + 8.2961*$t*$t*$t - 265.5*$t*$t + 3739*$t - 19639} elsif ($t<31) {$r=224} else {$r=0.8333*$t*$t*$t*$t - 108.76*$t*$t*$t + 5312.7*$t*$t -115145*$t + 934640}
		if ($t<3) {$g=8} elsif ($t<13) {$g=-0.7331*$t*$t + 35.011*$t - 93.933} elsif ($t<26) {$g=-0.076*$t*$t*$t*$t + 6.04567*$t*$t*$t - 175.04*$t*$t + 2178.2*$t - 9617.5} else {$g=-0.0758*$t*$t*$t*$t + 10.04*$t*$t*$t - 492.3*$t*$t + 10570*$t - 83619}
		if ($t<5) {$b=1.0648*$t*$t*$t - 10.04*$t*$t - 3.4299*$t + 238.34} elsif ($t<13) {$b = -0.1126*$t*$t + 18.593*$t + 12.978} elsif ($t<22) {$b=0.037*$t*$t*$t*$t - 3.1606*$t*$t*$t + 98.528*$t*$t -1353.1*$t + 7069.9} else {$b=8}
		$r = int($r); $g = int($g); $b = int($b);
		$out = [$r,$g,$b];
	} else {die "\n#### FATAL ERROR ####\ngradient type is not recognized: $gradient_type\n";}
	return $out;
}
