#!/usr/bin/perl -w
use strict;

die <<END unless @ARGV >= 3;

 *******************************************************************************************************************************
  Author: wangsc\@ksu.edu
  This script will calcualte the haplotype share statistics (PHS) as reported in the paper:																				
  A Nonparametric test reveals selection for rapid flowering in the arabidopsis genome, Christopher T et al, PLoS Biology 2006.
 
  Jun.29.2013 by SW:
  Modified from calcualte_PHS_v3.pl, to calcualte the mean and sd of the length of shared haplotype from all chromosmes.
  
  Bug fixed on 08.08.2013, SW
  
  Input files needed:
  File_1. Ordered mapping distance of SNPs for all chromosomes.
     Format:
     chr_id SNP_name genetic_distance
     1A wsnp_Ex_c12345 5.8
     1A wsnp_Ex_c34590 6.5
     ...
     1B wsnp_Ex_c3450 16.5
   
  
  File_2. Genotypping file for each accession. 
  ** Only genotyping data for SNPs presented in File_1 will be used for calculation.
     Format:
     accession_id wsnp_Ex_c12345 wsnp_Ex_c34590 .  .  ....
     acc1          A             B              B  NA ...
     acc2          A             B              A  A  ...
     ...
  
  Output files:
  <chromosome>_phs.txt
    
 
  Usage:
    perl $0
    <Mapping distance file>
    <Genotyping file>
    <chromosome>

  Example: perl $0  snps.dist  genotyping_all_snps.txt 1B
 
 *******************************************************************************************************************************
END

my ($mapped_marker_file, $genotype_file, $chromosome) = @ARGV ;

print_time("Reading marker file $mapped_marker_file");
my ($mapped_markers_ref, $marker_genetic_dist_href, $marker_chr_href) = read_marker_file($mapped_marker_file);

print_time("Reading genotype file $genotype_file");
my %genotype_marker = read_genotype_file($genotype_file, $mapped_markers_ref);

my %allele_freq_each_marker = count_allele_freq(\%genotype_marker);
my @accessions = keys %genotype_marker;
print STDERR "No. accessions: ", scalar @accessions, "\n";

my %haplotype_extent;
my %marker_blocks;
my @pairs;
my %block_length_hash;
my $count = 0;

print_time("Calculate pairwise haplotype extent ...");
my @markers = @$mapped_markers_ref;	
foreach my $acc_a (0..($#accessions-1))
{
	foreach my $acc_b ($acc_a+1 .. $#accessions)
	{
		if ($accessions[$acc_a]=~/:/ or $accessions[$acc_b]=~/:/){print STDERR "!! Accession names should NOT have \":\"\n"; exit}
		$count++;
		print_time("\tFinished $count pairs ..") unless $count%100;
		my $pair = join(":", (sort{$a cmp $b} ($accessions[$acc_a], $accessions[$acc_b])));
		push @pairs, $pair;
    my @block_length = ();
		
		&pairwise_haplotype_extent_calc
		(\%genotype_marker, \%haplotype_extent, 
		 \%marker_blocks, $chromosome,
		 \@markers, $marker_genetic_dist_href, 
		 $accessions[$acc_a], $accessions[$acc_b], 
		 $pair, $marker_chr_href, \@block_length);
        
    $block_length_hash{$pair} = [@block_length];
	}
}

# calculate genome_wide_Mean_Sigma of blocks
print_time("Calculate genome_wide_Mean_Sigma of blocks ...");
#my %mean_sd_hash = calculate_genome_wide_Mean_Sigma(\%haplotype_extent, \@pairs, $marker_genetic_dist_href, $mapped_markers_ref, $marker_chr_href);
my %mean_sd_hash = map {$_, [average($block_length_hash{$_}), stdev($block_length_hash{$_})] } @pairs;

# Calculate PHS for the selected chromosome
foreach my $chr(sort{$a cmp $b} unique(values %$marker_chr_href))
{
	my %PHS;
	my %average_block;
	next unless $chr eq $chromosome; # selected chromosome only;
	# get markers on the chromsome
	my @markers_chr;
	map{push @markers_chr, $_ if $marker_chr_href->{$_} eq $chr}@markers;
	
	print_time("Calculating PHS for chr $chr ...");		
	foreach my $index (0 .. $#markers_chr)
	{
		my $marker_name = $markers_chr[$index];
		print STDERR $index, "\t", $marker_name, "\n" if $index == 0;	  
		# calculate N Z(xjx)
		my ($N_sumZ_avg, $A_sumZ_avg, $B_sumZ_avg, $A_avg_block, $B_avg_block) = 
        sum_Z($marker_name, \%genotype_marker, \@pairs, \%marker_blocks, 
		          $index, $marker_genetic_dist_href, $mapped_markers_ref, \%mean_sd_hash, $chr);
		
		$PHS{$index} = [$A_sumZ_avg - $N_sumZ_avg, $B_sumZ_avg - $N_sumZ_avg];
		$average_block{$index} = [$A_avg_block, $B_avg_block];
	}
	
	# Output
	my $phs_out_file = "chromosome_".$chr . "_phs.txt";
	open (PHS, ">$phs_out_file") or die $!;
	print PHS join("\t", qw(Chr SNP_index SNP_name Genetic_dist PHS_A PHS_B Freq_A Freq_B Block_A Block_B)),"\n";
	
	foreach my $index (0 .. $#markers_chr)
	{
		next unless exists $PHS{$index};
		my $marker = $markers_chr[$index];
		my $genetic_dist = $marker_genetic_dist_href->{$marker};
		my @phs = @{$PHS{$index}};
		my @avg_left_right_block = @{$average_block{$index}};
		next if (not exists $allele_freq_each_marker{$marker});
		my @allele_freq = @{$allele_freq_each_marker{$marker}};
		print PHS join("\t",($chr, $index, $marker, $genetic_dist, @phs, @allele_freq, @avg_left_right_block)),"\n";
	}
	
	close PHS;
	print_time("Finished PHS calculation for chromosome $chr ...");
	print_time("Generate PHS resut file $phs_out_file for chromosome $chr");
}
print_time("All Done...");

# Subroutines

sub sum_Z
{
	my ($marker_name, $genotype_ref, $pairs_ref, 
	    $marker_blocks_ref, $index, $marker_dist_ref, $mapped_markers_ref, $mean_sd_hashref, $chr) = @_;
	
	my @markers_chr; #SW 08.06.2013
	map{push @markers_chr, $_ if $marker_chr_href->{$_} eq $chr}@$mapped_markers_ref; #SW 08.06.2013
	
	my $no_pairs = scalar @$pairs_ref;
	
	my @sum_z = (0,0,0);
	my @avg_left_right_block_A = (0,0);
	my @avg_left_right_block_B = (0,0);
	my @count = (0,0,0);
	foreach my $pair (@$pairs_ref)
	{
		next unless exists $marker_blocks_ref->{$chr}{$index}{$pair};
		next unless exists $mean_sd_hashref->{$pair};
		my ($mean, $sd) = @{$mean_sd_hashref->{$pair}};
		unless (defined $sd){print STDERR $pair, "\n"}
		
		my ($acc1, $acc2) = split /:/, $pair;
		my $genotype1 = $genotype_ref->{$acc1}{$marker_name};
		my $genotype2 = $genotype_ref->{$acc2}{$marker_name};
		next if $genotype1 eq "NA" or $genotype2 eq "NA";
		if($genotype1 ne $genotype2){print STDERR "Err, strange ...\n"; next}
	
		my $block = $marker_blocks_ref->{$chr}{$index}->{$pair};
		next unless defined $block;
		my ($start, $end) = (split/_/, $block);
		my ($left_dist, $right_dist) = sort{$a<=>$b} ($marker_dist_ref->{$markers_chr[$start]}, $marker_dist_ref->{$markers_chr[$end]}); #SW 08.06.2013
		my $dist = $right_dist - $left_dist;
		$sum_z[0] += $sd>0?($dist - $mean)/$sd:0;
		$count[0]++;
		if($genotype1 eq "A"){$count[1]++; $sum_z[1] += $sd>0?($dist - $mean)/$sd:0; $avg_left_right_block_A[0] += $left_dist; $avg_left_right_block_A[1] += $right_dist}
		if($genotype1 eq "B"){$count[2]++; $sum_z[2] += $sd>0?($dist - $mean)/$sd:0; $avg_left_right_block_B[0] += $left_dist; $avg_left_right_block_B[1] += $right_dist}
		
		if($marker_name eq "4385276_6as:2126")
		{
			print STDERR "genotype1: ", $genotype1, "\n", "genotype2: ", $genotype2, "\n",
			             "sum_z: ", join(" ", @sum_z),"\n", 
			             "avg_A: ", join(" ", @avg_left_right_block_A), "\n", 
			             "avg_B: ", join(" ", @avg_left_right_block_B),"\n" ;
		}
	}
	
	map{$sum_z[$_] =  $sum_z[$_]/$count[$_] if $count[$_]>0}0..$#sum_z;
	map{$avg_left_right_block_A[$_] = $avg_left_right_block_A[$_]/$count[1] if $count[1]>0} 0 .. $#avg_left_right_block_A;
	map{$avg_left_right_block_B[$_] = $avg_left_right_block_B[$_]/$count[2] if $count[2]>0} 0 .. $#avg_left_right_block_B;
	
	return (@sum_z, join("_", @avg_left_right_block_A), join("_", @avg_left_right_block_B));
}


sub CN2
{
	my $N = shift;
	# calculate C(N,2);
	# C(4,2) = 6; C(5,2) = 10;
	return ($N-1)*$N/2;
}


sub calculate_genome_wide_Mean_Sigma
# return (mean, stand_deviation);
{
	my ($hap_extend_freq, $accession_pairs_ref, $marker_genetic_dist_href, $mapped_markers_ref, $marker_chr) = @_;
	my @acc_pairs = @$accession_pairs_ref;
	my %return;

	my @markers = @{$mapped_markers_ref};
	foreach my $pair (@acc_pairs)
	{
		my @block_length;
		foreach my $chr (unique(values %$marker_chr))
		{
			next unless defined $hap_extend_freq->{$chr}{$pair};
			my @blocks = @{$hap_extend_freq->{$chr}{$pair}};
			foreach (@blocks)
			{
				my @two_markers = @markers[(split/_/,$_)];
				my $dist  = abs($marker_genetic_dist_href->{$two_markers[0]} - $marker_genetic_dist_href->{$two_markers[1]});
				push @block_length, $dist;
			}			
		}

		$return{$pair} = [average(\@block_length), stdev(\@block_length)];
	}
	
	return %return;
}

sub pairs
{
	my $arr_ref = shift;
	return unless @$arr_ref > 1;
	my @return;
	foreach my $index (0..(scalar @$arr_ref)-2)
	{
		foreach ($index .. (scalar @$arr_ref)-1)
		{
			push @return, join(":", sort {$a cmp $b} ($arr_ref->[$index], $arr_ref->[$_]));
		}
	}
	return @return;
}


sub read_marker_file
{
	my $file = shift;
	my @markers_each_chr;
	my %distance_each_marker;
	my %chromosome;
	open (IN, $file) or die;
	my $line = 0;
	while(<IN>)
	{
		chomp;
		$line++;
		next if $line==1; #/chr_id/i; # skip first line
		my ($chr, $marker, $distance) = split /\s+/, $_;
		if(exists $chromosome{$marker}){print STDERR  "\t!!! Marker $marker presents in multiple chromosomes. !!!\n"}
    else{$chromosome{$marker} = $chr;}
		push @markers_each_chr, $marker;
		$distance_each_marker{$marker} = $distance;
	}
	close IN;
	return(\@markers_each_chr, \%distance_each_marker, \%chromosome);
}

sub read_genotype_file
{
	my ($genotype_file, $markers_each_chr) = @_;
	my %return;
	my @mapped_markers = @$markers_each_chr;
	my %mapped_markers_hash = map{$_, 1} @mapped_markers;
	open (IN, $genotype_file) or die;
	my $line = 0;
	my @total_markers;
	my @index_for_mapped_markers;
	while(<IN>)
	{
		next if /^\s+$/;
		$line++;
		chomp;
		s/\s+$//;
		my @t = split /\s+/,$_;
		if ($line == 1) # get marker names from the first line
		{
			@total_markers = @t;
			foreach (1..$#t)
			{
				push @index_for_mapped_markers, $_ if exists $mapped_markers_hash{$t[$_]}
			}
			next;
		}
		my $acc = $t[0];
		foreach (@index_for_mapped_markers)
		{
			die $_, "\t", scalar @total_markers, "\n" unless defined $total_markers[$_];
			$return{$acc}{$total_markers[$_]} = $t[$_];
		}		
	}
	close IN;
	return %return;
}

sub count_allele_freq
{
	my $genotype_hashref = shift;
	my %return;
	my $num_accessions = scalar (keys %$genotype_hashref);
	foreach my $acc (keys %$genotype_hashref)
	{
		foreach my $mrk (keys %{$genotype_hashref->{$acc}})
		{
			my $allele = $genotype_hashref->{$acc}{$mrk};
			$return{$mrk}{$allele}++;
		}
	}
	
	foreach my $mrk (keys %return)
	{
		my @alleles = keys %{$return{$mrk}};
		@alleles = sort {$a cmp $b} @alleles;
		$return{$mrk}{"A"} = 0 unless exists $return{$mrk}{"A"};
		$return{$mrk}{"B"} = 0 unless exists $return{$mrk}{"B"};
		my $total = $return{$mrk}{"A"} + $return{$mrk}{"B"};
		if ($total > 0)
		{
			$return{$mrk} = [$return{$mrk}{"A"}/$total, $return{$mrk}{"B"}/$total];
		}
		else
		{
			$return{$mrk} = [0, 0];
		}
		
	}
	return %return;
}

sub pairwise_haplotype_extent_calc
{
	my ($gn_marker, $hap_ext_freq, $marker_blocks, $chr_tobe_calculated, $markers, 
        $marker_genetic_dist_href, $acc1, $acc2, $pair, $marker_chr, $block_length_ref) = @_;
	my @chrs = sort{$a cmp $b} unique(values %$marker_chr);
	
	foreach my $chr(@chrs)
	{
		my @status_record;
		my @markers_chr;
		map{push @markers_chr, $_ if $marker_chr_href->{$_} eq $chr}@markers;
		foreach my $index (0 ..  $#markers_chr)
		{
			my $mk = $markers_chr[$index];
      #print STDERR "pairwise: ", $index, "\t", $mk, "\n" if $index == 0;
			$gn_marker->{$acc1}->{$mk} = "NA" unless exists $gn_marker->{$acc1}->{$mk};
			$gn_marker->{$acc2}->{$mk} = "NA" unless exists $gn_marker->{$acc2}->{$mk};
			if ($gn_marker->{$acc1}->{$mk} eq $gn_marker->{$acc2}->{$mk} and $gn_marker->{$acc2}->{$mk} ne "NA")
			{
				$status_record[$index] = 1;
			}
			elsif($gn_marker->{$acc1}->{$mk} eq "NA" or $gn_marker->{$acc2}->{$mk} eq "NA")
			{
				$status_record[$index] = -1;
			}
			else
			{
				$status_record[$index] = 0;
			}
			die "status_record $index not defiend\n" unless defined $status_record[$index];
		}
		
		my @blocks = calculate_blocks(@status_record);
    map{
            my @p=split /_/,$_; 
            push @{$block_length_ref}, abs($marker_genetic_dist_href->{$markers_chr[$p[0]]} - $marker_genetic_dist_href->{$markers_chr[$p[1]]});
    }@blocks;
        
    if($chr_tobe_calculated eq $chr)
    {
	    	$hap_ext_freq->{$chr}->{$pair} = [@blocks] if @blocks > 0;
	    	
	    	foreach my $blk (@blocks)
	      {
	    		my ($start, $end) = split /_/, $blk;
	    		foreach ($start .. $end)
	    		{
	    			$marker_blocks->{$chr}->{$_}{$pair} = $blk;
	    		}
	    	}
     }
	}

}

sub calculate_blocks
{
	my @records = @_;
	my @return;
	my $pre; 
	foreach my $index (0..$#records)
	{
		if($records[$index] == 1)
		{
			$pre = $index unless defined $pre;
			if($index == $#records)
			{
				push @return, [$pre, $index] unless $pre == $index;
			}
		}
		elsif($records[$index] == 0)
		{
			if(defined $pre)
			{
				push @return, [$pre, $index-1] unless $pre == $index-1;
				undef $pre;
			}
		}
		elsif($records[$index] == -1)
		{
			if(defined $pre)
			{
				if(($index < $#records and $records[$index+1] != 1) or ($index == $#records))
				{
					push @return, [$pre, $index-1] unless $pre == $index-1;
					undef $pre;
				}
			}
		}
	}
	@return = map{join("_", @$_)}@return;
	#print join("", @records), "\n";
	return @return;
}

sub stdev {
	my $arr_ref = $_[0];

	my $n = 0;
	my $sumsq = 0;
	my $sum = 0;
	foreach (@$arr_ref) {
		next unless (defined($_) && ($_ !~ /^$/));
		$n++;
		$sumsq += $_ * $_;
		$sum += $_;
	}
	if ($n <= 1) { return 0; }
	my $stdev = (($n * $sumsq) - ($sum * $sum))/($n * ($n - 1));
	$stdev = ($stdev < 0) ? 0 : $stdev;

	return (sqrt($stdev));
}

sub average {
	my $count = &count($_[0]);
	if ($count == 0) {
		return 0;
	} else {
		return (&sum($_[0]) / $count);
	}
}

sub sum {
	my $arr_ref = $_[0];

	my $sum = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$sum += $_;
	}
	return $sum;
}

sub count {
	my $arr_ref = $_[0];
	
	my $num = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$num++;
	}
	return $num;
}

sub unique
{
	my %temp = map{$_, 1}@_;
	return keys %temp;
}

sub print_time
{
	my $str = shift;
	my $time = localtime(time);
	print STDERR join("  ", ($time, $str)), "\n";
}
