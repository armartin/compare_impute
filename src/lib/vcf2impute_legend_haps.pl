#!/usr/bin/perl
#
# [ vcf2impute_legend_haps.pl ]
#
# This script takes a VCF file and converts it into a haplotype+legend file
# pair in IMPUTE format.
#
# *** Copyright Bryan Howie, 2013 ***
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

sub InitializeSNPAllelesTable();
sub InitializeDtypeLookupTable();
sub ReadSampleTableFile($\%\%\%);
sub PrintLegendHapsFiles($$);
sub ParseCommandLine();

ParseCommandLine();

# command line variables
my %arg;
my $vcf_file = $arg{-vcf};                # VCF file to process; should contain phased haplotypes
my $new_basename = $arg{-leghap};         # basename of new files to print; suffixes and gzipping are automatic
my $chrom = $arg{-chr};                   # chromosome to include in output files
my $sample_table_file = $arg{-samp_tab};  # file containing 3 columns: sample ID, pop ID, group ID
my $keep_group = $arg{-keep_group};       # string of a single group to include in the output
my $var_counts_file = $arg{-var_counts};  # file to print with counts of each variant type (e.g., 'SNP')
my $min_maf = $arg{-min_maf};             # smallest allowed MAF; applied only to -keep_group when active
my $start_pos = $arg{-start};             # first position to include in output files
my $end_pos = $arg{-end};                 # last position to include in output files
my $no_mono = $arg{-no_mono};             # flag: omit monomorphic sites from processed files
my $snps_only = $arg{-snps_only};         # flag: omit any sites that are not biallelic SNPs
my $is_chrX_nonpar = $arg{-chrX_nonpar};  # flag: we are analyzing a non-PAR region of chrX
my $keep_annot = $arg{-keep_annot};       # flag: keep certain annotations from VCF (e.g. VT, SNPSOURCE)
my $old_var_ids = $arg{-old_var_ids};     # flag: use old [chr]-[pos] convention for variants w/o IDs

# global variables
my $line;
my (%snp_alleles, %phased_gt_to_probs_diploid, %phased_gt_to_probs_haploid,
    %sample_to_ancestry_group, %n_haps_per_group, %sample_ploidy);
my $chrX_par1_start_hg19 = 60001;
my $chrX_par1_end_hg19 = 2699520;
my $chrX_par2_start_hg19 = 154931044;
my $chrX_par2_end_hg19 = 155260560;

# file paths and temporary file names


MAIN: {

  #
  InitializeSNPAllelesTable();

  #
  InitializeDtypeLookupTable();

  #
  ReadSampleTableFile($sample_table_file, %sample_to_ancestry_group, %n_haps_per_group, %sample_ploidy) if ($sample_table_file ne "");

  #
  PrintLegendHapsFiles($vcf_file, $new_basename);

}



###################################################################
#                            FUNCTIONS                            #
###################################################################

#
sub InitializeSNPAllelesTable() {
  $snp_alleles{'A'} = 1;
  $snp_alleles{'a'} = 1;
  $snp_alleles{'C'} = 1;
  $snp_alleles{'c'} = 1;
  $snp_alleles{'G'} = 1;
  $snp_alleles{'g'} = 1;
  $snp_alleles{'T'} = 1;
  $snp_alleles{'t'} = 1;
}


#
sub InitializeDtypeLookupTable() {
  $phased_gt_to_probs_diploid{'.|.'} = '? ?';
  $phased_gt_to_probs_diploid{'0|0'} = '0 0';
  $phased_gt_to_probs_diploid{'0|1'} = '0 1';
  $phased_gt_to_probs_diploid{'1|0'} = '1 0';
  $phased_gt_to_probs_diploid{'1|1'} = '1 1';

  $phased_gt_to_probs_haploid{'.|.'} = '? -';
  $phased_gt_to_probs_haploid{'0|0'} = '0 -';
  $phased_gt_to_probs_haploid{'0|1'} = '? -';
  $phased_gt_to_probs_haploid{'1|0'} = '? -';
  $phased_gt_to_probs_haploid{'1|1'} = '1 -';
  $phased_gt_to_probs_haploid{'.'}   = '? -';
  $phased_gt_to_probs_haploid{'0'}   = '0 -';
  $phased_gt_to_probs_haploid{'1'}   = '1 -';
}


# read sample table; first column gives sample ID, second column gives
# population ID, and third column gives group ID (e.g., continent)
sub ReadSampleTableFile($\%\%\%) {
  my ($samp_tab_file, $sample_to_group_table, $group_hap_counts_table, $ploidy_table) = @_;
  my @line_info;
  my $ploidy;

  open(SAMP_TAB, $samp_tab_file) || die $!;

  while (defined ($line = <SAMP_TAB>)) {
    @line_info = split /\s+/, $line;

    if (scalar(@line_info) > 3 && $chrom eq "X" && $is_chrX_nonpar) {
      if ($line_info[3] == 1) {
	$ploidy = 1;  # male
      }
      else {
	$ploidy = 2;  # female
      }
    }
    else {
      $ploidy = 2;  # default
    }
    $$ploidy_table{$line_info[0]} = $ploidy;

    $$sample_to_group_table{$line_info[0]} = $line_info[2];
    if (!exists($$group_hap_counts_table{$line_info[2]})) {
      $$group_hap_counts_table{$line_info[2]} = $ploidy;
    }
    else {
      $$group_hap_counts_table{$line_info[2]} += $ploidy;
    }
  }

  close(SAMP_TAB);
}


#
sub PrintLegendHapsFiles($$)
{
  my ($infile, $basename) = @_;
  my $leg_file = $basename.".legend.gz";
  my $hap_file = $basename.".hap.gz";
  my $samp_file = $basename.".sample_list";
  my (@line_info, @pop_groups, %nonref_allele_counts_by_group, %var_counts);

  # open input VCF file, allowing for possible gzipping
  my $is_in_gz = ($infile =~ m/.gz$/);
  open(VCF, ($is_in_gz ? "gzip -dc $infile |" : $infile)) || die "Can't open $infile: $!\n";

  # create a handle for writing a new gzipped .leg file
  open(LEG, " | gzip -c > $leg_file") || die "Can't open $leg_file: $!\n";

  # create a handle for writing a new gzipped .haps file
  open(HAP, " | gzip -c > $hap_file") || die "Can't open $hap_file: $!\n";

  # find last header line and use it to print a sample list
  my $pos_tag_idx = -1;
  my $id_tag_idx = -1;
  my $ref_tag_idx = -1;
  my $alt_tag_idx = -1;
  my $filter_tag_idx = -1;
  my $info_tag_idx = -1;
  my $format_tag_idx = -1;
  my @sample_list;

  while (defined ($line = <VCF>)) {
    chomp $line;
    @line_info = split /\t/, $line;
    if (scalar(@line_info) > 0 && $line_info[0] eq "#CHROM") {
      for my $i (1..$#line_info) {
	if ($line_info[$i] eq "POS") {
	  $pos_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "ID") {
	  $id_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "REF") {
	  $ref_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "ALT") {
	  $alt_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "FILTER") {
	  $filter_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "INFO") {
	  $info_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "FORMAT") {
	  $format_tag_idx = $i;
	  last;
	}
      }

      # no POS tag found; something must be wrong
      if ($pos_tag_idx == -1) {
	die "\nERROR: No POS tag found.\n\n";
      }

      # no ID tag found; something must be wrong
      if ($id_tag_idx == -1) {
	die "\nERROR: No ID tag found.\n\n";
      }

      # no REF tag found; something must be wrong
      if ($ref_tag_idx == -1) {
	die "\nERROR: No REF tag found.\n\n";
      }

      # no ALT tag found; something must be wrong
      if ($alt_tag_idx == -1) {
	die "\nERROR: No ALT tag found.\n\n";
      }

      # no FILTER tag found; something must be wrong
      if ($filter_tag_idx == -1) {
	die "\nERROR: No FILTER tag found.\n\n";
      }

      # no FORMAT tag found; something must be wrong
      if ($format_tag_idx == -1) {
	die "\nERROR: No FORMAT tag found.\n\n";
      }

      # if we reach this point, we must have successfully found the FORMAT tag;
      # all subsequent entries should be sample IDs
      open(SAMP, "> $samp_file") || die $!;
      for my $j (($format_tag_idx+1)..$#line_info) {
	push @sample_list, $line_info[$j];
	print SAMP "$line_info[$j]\n" if ($keep_group eq "" || $keep_group eq $sample_to_ancestry_group{$line_info[$j]});
      }
      close(SAMP);

      last; # done reading header
    }
  }

  # print legend header
  print LEG "id position a0 a1";
  if ($keep_annot) {
    print LEG " type source rsq";
  }
  if (scalar(keys(%n_haps_per_group)) > 0) {
    @pop_groups = sort(keys(%n_haps_per_group));
    foreach my $group (@pop_groups) {
      print LEG " ".lc($group).".aaf" if ($keep_group eq "" || $keep_group eq $group);
    }
    foreach my $group (@pop_groups) {
      print LEG " ".lc($group).".maf" if ($keep_group eq "" || $keep_group eq $group);
    }
  }
  print LEG "\n";

  # process SNP information
  my $n_snps = 0;
  my $n_filtered_snps = 0;
  my $n_multiallelic_snps = 0;
  my $n_monomorphic_snps = 0;
  my $n_non_snp_variants = 0;
  my $n_low_maf_snps = 0;
  my $is_chrom = 0;  # have we reached the chromosome of interest yet?
  my $var_id = "";
  my $var_type = ".";  # options are SNP, INDEL, and SV
  my $var_source = ".";  # options are LOWCOV and EXOME
  my $rsq = ".";  # genotype imputation quality from MaCH/Thunder
  my (@dtype_info, @alleles);
  while (defined ($line = <VCF>)) {
    chomp $line;
    @line_info = split /\t/, $line;

    if (++$n_snps % 5000 == 0) {
      print STDERR "chr$line_info[0] -- $line_info[$pos_tag_idx]\n";
    }

    $is_chrom = 1 if ($line_info[0] eq $chrom);
    if ($line_info[0] ne $chrom) {
      last if ($is_chrom);  # assumes that all VCF rows for this chrom are contiguous
      next;
    }

    if ($line_info[$filter_tag_idx] ne "PASS" && $line_info[$filter_tag_idx] ne ".") {
      ++$n_filtered_snps;
      next;
    }

    if ($line_info[$alt_tag_idx] =~ ',') {
      ++$n_multiallelic_snps;
      next;
    }

    if ($line_info[$alt_tag_idx] eq '.' && $no_mono) {
      ++$n_monomorphic_snps;
      next;
    }

    if ($snps_only &&
	!(exists($snp_alleles{$line_info[$ref_tag_idx]}) && exists($snp_alleles{$line_info[$alt_tag_idx]})) ) {
      ++$n_non_snp_variants;
      next;
    }

    next if ($line_info[$pos_tag_idx] < $start_pos);
    last if ($line_info[$pos_tag_idx] > $end_pos);

    # special handling of the non-pseudoautosomal region of chromosome X
    next if ($is_chrX_nonpar && $line_info[$pos_tag_idx] <= $chrX_par1_end_hg19);
    last if ($is_chrX_nonpar && $line_info[$pos_tag_idx] >= $chrX_par2_start_hg19);

    # if we are computing allele freqs within ancestral groups, initialize the allele counts here
    if (scalar(@pop_groups) > 0) {
      foreach my $group (@pop_groups) {
	$nonref_allele_counts_by_group{$group} = 0;
      }
    }

    if ($keep_annot || $var_counts_file ne "") {
      my @annot = split /;/, $line_info[$info_tag_idx];
      $var_type = '.';
      $var_source = '.';
      $rsq = '.';
      foreach (@annot) {
	my @curr_annot = split /=/, $_;
	$var_type = $curr_annot[1] if ($curr_annot[0] eq "VT");
	$var_source = $curr_annot[1] if ($curr_annot[0] eq "SNPSOURCE");
	$rsq = $curr_annot[1] if ($curr_annot[0] eq "RSQ");
      }
      $var_source = "LOWCOV" if ($var_source eq '.');  # as of 1000G Phase 1, the variant source is not identified for INDEL/SV, which all come from low-cov data
    }

    # parse haplotype alleles; also store group-level allele counts if computing allele freqs
    my $is_first_samp = 1;  # boolean
    my $j = 0;
    my $hap_alleles_str = "";
    for my $i (($format_tag_idx+1)..$#line_info) {
      my $curr_group = "";
      $curr_group = $sample_to_ancestry_group{$sample_list[$j]} if (scalar(@pop_groups) > 0);
      if ($curr_group ne "" && $keep_group ne "" && $keep_group ne $curr_group) {
	++$j;
	next;
      }

      @dtype_info = split /\:/, $line_info[$i];
      if ($is_chrX_nonpar && $sample_ploidy{$sample_list[$j]} == 1) {  # male on chrX nonpar; haploid
	if (exists($phased_gt_to_probs_haploid{$dtype_info[0]})) {
	  if ($dtype_info[0] eq "0|1" || $dtype_info[0] eq "1|0") {
	    die "\nERROR: Diplotype |$dtype_info[0]| (i=$i) at site $line_info[$pos_tag_idx] is not in accepted format; hemizygous males cannot have heterozygous genotypes.\n\n";
	  }
	  $hap_alleles_str .= ($is_first_samp ? "" : " ")."$phased_gt_to_probs_haploid{$dtype_info[0]}";
	}
	else {
	  die "\nERROR: Diplotype |$dtype_info[0]| (i=$i) at site $line_info[$pos_tag_idx] is not in accepted format.\n\n";
	}
      }
      else {  # dizygous individual
	if (exists($phased_gt_to_probs_diploid{$dtype_info[0]})) {
	  $hap_alleles_str .= ($is_first_samp ? "" : " ")."$phased_gt_to_probs_diploid{$dtype_info[0]}";
	}
	else {
	  die "\nERROR: Diplotype |$dtype_info[0]| (i=$i) at site $line_info[$pos_tag_idx] is not in accepted format.\n\n";
	}
      }
      $is_first_samp = 0;

      # if computing allele freqs, update allele count here
      if ($curr_group ne "") {
	if ($is_chrX_nonpar && $sample_ploidy{$sample_list[$j]} == 1) {  # male on chrX nonpar; haploid
	  if ($dtype_info[0] eq '1|1') {
	    $nonref_allele_counts_by_group{$curr_group} += 1;
	  }
	}
	else {  # diploid individual
	  if ($dtype_info[0] eq '0|1' || $dtype_info[0] eq '1|0') {
	    ++$nonref_allele_counts_by_group{$curr_group};
	  }
	  elsif ($dtype_info[0] eq '1|1') {
	    $nonref_allele_counts_by_group{$curr_group} += 2;
	  }
	}
      }
      ++$j;
    }

    # create a string to store legend info beyond the foundational four columns
    my $legend_str = "";

    # if keeping annotations from the VCF, add the info from this variant to the legend string
    if ($keep_annot) {
      $legend_str .= " $var_type $var_source $rsq";
    }

    # if computing allele freqs, do so here based on the allele counts tallied
    # above; also store the corresponding strings for printing in the legend file below
    my $curr_max_maf = 0;  # largest MAF across groups at this variant; used for -min_maf filtering
    if (scalar(@pop_groups) > 0) {
      ## 'alternate' allele freqs
      foreach my $group (@pop_groups) {
	next if ($keep_group ne "" && $keep_group ne $group);
	$legend_str .= " ".sprintf("%.4f", $nonref_allele_counts_by_group{$group} / $n_haps_per_group{$group});
      }

      ## minor allele freqs
      foreach my $group (@pop_groups) {
	next if ($keep_group ne "" && $keep_group ne $group);
	my @allele_counts;
	push @allele_counts, $nonref_allele_counts_by_group{$group};
	push @allele_counts, $n_haps_per_group{$group} - $nonref_allele_counts_by_group{$group};
	@allele_counts = sort {$a <=> $b} @allele_counts;
	my $n_minor_alleles = $allele_counts[0];
	my $curr_maf = $n_minor_alleles / $n_haps_per_group{$group};
	$legend_str .= " ".sprintf("%.4f", $curr_maf);
	$curr_max_maf = $curr_maf if ($curr_maf > $curr_max_maf);
      }
    }

    # if the highest MAF of this variant in any population group (or just the -keep_group
    # when that option is active) is below the -min_maf, omit it from the output
    if ($curr_max_maf < $min_maf) {
      ++$n_low_maf_snps;
      next;
    }

    #
    if ($var_counts_file ne "") {
      if (exists($var_counts{$var_type})) {
	++$var_counts{$var_type};
      }
      else {
	$var_counts{$var_type} = 1;
      }
    }

    # if we reach this point of the chain of execution, the variant has passed all of
    # the prescribed filters, so we can print its legend and haplotype file entries
    $var_id = $line_info[$id_tag_idx];
    my $ref_allele = $line_info[$ref_tag_idx];
    my $alt_allele = $line_info[$alt_tag_idx];
    $ref_allele = "-" if ($ref_allele eq "<DEL>");
    $alt_allele = "-" if ($alt_allele eq "<DEL>");
    if ($var_id eq '.') {
      if ($old_var_ids) {
	$var_id = "$chrom"."-$line_info[$pos_tag_idx]";
      }
      else {
	$var_id = "chr$chrom".":$line_info[$pos_tag_idx]";
	if (exists($snp_alleles{$ref_allele}) && length($alt_allele) > 1) {
	  $var_id .= ":I";
	}
	elsif (length($ref_allele) > 1 && exists($snp_alleles{$alt_allele})) {
	  $var_id .= ":D";
	}
	elsif ($ref_allele eq '-') {
	  $var_id .= ":I";
	}
	elsif ($alt_allele eq '-') {
	  $var_id .= ":D";
	}
      }
    }
    print LEG "$var_id $line_info[$pos_tag_idx] $ref_allele $alt_allele".$legend_str."\n";
    print HAP "$hap_alleles_str\n";
  }

  # print number of SNPs removed by FILTER tag
  print STDERR "\n--$n_filtered_snps SNPs were removed by the FILTER tag.\n";
  print STDERR "\n--$n_multiallelic_snps SNPs were removed for having more than one allele in the ALT column.\n";
  print STDERR "\n--$n_monomorphic_snps SNPs were removed for having no defined allele in the ALT column.\n" if ($no_mono);
  print STDERR "\n--$n_non_snp_variants variants were removed for having non-SNP alleles.\n" if ($snps_only);
  print STDERR "\n--$n_low_maf_snps variants were removed for having MAFs below the -min_maf threshold.\n" if ($n_low_maf_snps > 0);
  print STDERR "\n";

  close(VCF);
  close(LEG);
  close(HAP);

  # if requested, print a file showing the count for each type of variant encountered
  if ($var_counts_file ne "") {
    open(COUNTS, " > $var_counts_file") || die "Can't open $var_counts_file: $!\n";
    my @sorted_var_types = sort(keys(%var_counts));
    foreach (@sorted_var_types) {
      print COUNTS "$_\t$var_counts{$_}\n";
    }
    close(COUNTS);
  }
}


# parse the command line
sub ParseCommandLine()
{
  my $usage = "\nvcf2impute_legend_haps.pl\n";
  $usage .= "  [-vcf]          VCF file to process; should contain phased haplotypes\n";
  $usage .= "  [-leghap]       basename of new files to print; suffixes and gzipping are auto\n";
  $usage .= "  [-chr]          chromosome to include in output files, in  (chr)[1-22,X]\n";
  $usage .= "  <-samp_tab>     file containing 3+ columns: sample ID, pop ID, group ID, ...\n";
  $usage .= "  <-keep_group>   string of a single group to include in the output\n";
  $usage .= "  <-var_counts>   file to print with counts of each variant type (e.g., 'SNP')\n";
  $usage .= "  <-min_maf>      smallest allowed MAF; applied only to -keep_group when active\n";
  $usage .= "  <-start>        first position to include in output files\n";
  $usage .= "  <-end>          last position to include in output files\n";
  $usage .= "  <-no_mono>      flag: omit monomorphic sites from processed files\n";
  $usage .= "  <-snps_only>    flag: omit any sites that are not biallelic SNPs\n";
  $usage .= "  <-chrX_nonpar>  flag: we are analyzing a non-PAR region of chrX\n";
  $usage .= "  <-keep_annot>   flag: keep certain annotations from VCF (e.g. VT, SNPSOURCE)\n";
  $usage .= "  <-old_var_ids>  flag: use old [chr]-[pos] convention for variants w/o IDs\n";
  $usage .= "\n";
  $usage .= "  (args in square brackets required; args in pointy brackets optional)\n";
  $usage .= "\n";

  # default values
  $arg{-vcf} = "";
  $arg{-leghap} = "";
  $arg{-chr} = "";
  $arg{-samp_tab} = "";
  $arg{-keep_group} = "";
  $arg{-var_counts} = "";
  $arg{-min_maf} = 0;
  $arg{-start} = 0;
  $arg{-end} = 1e9;
  $arg{-no_mono} = 0;  # boolean variable; false unless flag provided
  $arg{-snps_only} = 0;  # boolean variable; false unless flag provided
  $arg{-chrX_nonpar} = 0;  # boolean variable; false unless flag provided
  $arg{-keep_annot} = 0;  # boolean variable; false unless flag provided
  $arg{-old_var_ids} = 0;  # boolean variable; false unless flag provided

  # parse the command line
  for (my $i = 0; $i <= $#ARGV; ++$i) {
    if ($ARGV[$i] =~ /^-/) {
      if ($ARGV[$i] eq "-no_mono" || $ARGV[$i] eq "-snps_only" || $ARGV[$i] eq "-chrX_nonpar" ||
	  $ARGV[$i] eq "-keep_annot" || $ARGV[$i] eq "-old_var_ids") {
	$arg{$ARGV[$i]} = 1;  # set boolean to 'true'
      }
      else {
	$arg{$ARGV[$i]} = $ARGV[$i+1];
      }
    }
  }

  die ($usage) if ($arg{-vcf} eq "");
  die ($usage) if ($arg{-leghap} eq "");
  die ($usage) if ($arg{-chr} eq "");

  # the -no_mono flag implies that we should set a -min_maf, so apply
  # this here if a -min_maf was not given on the command line
  if ($arg{-no_mono} && $arg{-min_maf} == 0) {
    $arg{-min_maf} = 0.0000000001;
  }
}
