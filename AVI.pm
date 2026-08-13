=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

AVI - annotate variants with AVI scores

=head1 SYNOPSIS

 mv AVI.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin AVI,file=/FULL_PATH/avi_scores.tsv.gz
 ./vep -i variations.vcf --plugin AVI,file=/FULL_PATH/avi_scores.tsv.gz,raw=1
 ./vep -i variations.vcf --plugin AVI,file=/FULL_PATH/avi_scores.tsv.gz,raw=1,quantile=1

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves AVI scores for single nucleotide
 variants from a tabix-indexed, bgzip-compressed TSV file.

 The input file must have the following tab-separated header and a matching
 .tbi index:

 #CHROM  POS  REF  ALT  ag_pathogenicity  phred  quantile

 The PHRED score is reported by default as AVI_PHRED. The optional raw=1
 parameter additionally reports AVI_AG_PATHOGENICITY, and quantile=1 reports
 AVI_QUANTILE.

 The plugin annotates single nucleotide variants only.

 AVI scores:

 ag_pathogenicity (approximately -1.3 to +4.0)
 Raw score. Higher = more deleterious.

 quantile (0.0 to 1.0)
 Tail probability: fraction of all genome-wide SNVs with an equal or higher
 raw score. Lower = more deleterious.

 phred (0 to approximately 80)
 PHRED-scaled quantile. A phred of 20 means the variant is among the top 1%
 most deleterious substitutions genome-wide. Higher = more deleterious.

=cut

package AVI;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my %SCORE_INFO = (
  phred => {
    name => 'AVI_PHRED',
    description => 'PHRED-scaled quantile (0 to approximately 80). A phred of 20 means the variant is among the top 1% most deleterious substitutions genome-wide. Higher = more deleterious.',
  },
  ag_pathogenicity => {
    name => 'AVI_AG_PATHOGENICITY',
    description => 'Raw score (approximately -1.3 to +4.0). Higher = more deleterious.',
  },
  quantile => {
    name => 'AVI_QUANTILE',
    description => 'Tail probability (0.0 to 1.0): fraction of all genome-wide SNVs with an equal or higher raw score. Lower = more deleterious.',
  },
);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  my ($positional_file) = grep { index($_, '=') == -1 } @{$self->params};
  my $file = $params->{file} || $positional_file;

  die "ERROR: No AVI file specified\n" unless $file;

  $self->{include_raw} = _enabled($params->{raw});
  $self->{include_quantile} = _enabled($params->{quantile});
  $self->add_file($file);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header = (
    $SCORE_INFO{phred}->{name} => $SCORE_INFO{phred}->{description},
  );

  $header{$SCORE_INFO{ag_pathogenicity}->{name}} = $SCORE_INFO{ag_pathogenicity}->{description}
    if $self->{include_raw};
  $header{$SCORE_INFO{quantile}->{name}} = $SCORE_INFO{quantile}->{description}
    if $self->{include_quantile};

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my ($start, $end) = ($vf->{start}, $vf->{end});
  my $ref = $vf->ref_allele_string;
  my $allele = $tva->variation_feature_seq;

  return {} unless defined($start) && defined($end) && $start == $end;
  return {} unless defined($ref) && $ref =~ /^[ACGT]$/;
  return {} unless defined($allele) && $allele =~ /^[ACGT]$/;

  foreach my $record (@{$self->get_data($vf->{chr}, $start, $end)}) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref,
        alts   => [$allele],
        pos    => $start,
        strand => $vf->strand,
      },
      {
        ref  => $record->{ref},
        alts => [$record->{alt}],
        pos  => $record->{start},
      }
    );

    return $record->{result} if @{$matches};
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $ref, $alt, $raw, $phred, $quantile) = split /\t/, $line;

  return unless defined($quantile);

  my %result = (
    $SCORE_INFO{phred}->{name} => $phred,
  );

  $result{$SCORE_INFO{ag_pathogenicity}->{name}} = $raw if $self->{include_raw};
  $result{$SCORE_INFO{quantile}->{name}} = $quantile if $self->{include_quantile};

  return {
    chr    => $chr,
    start  => $start,
    end    => $start,
    ref    => $ref,
    alt    => $alt,
    result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub _enabled {
  my $value = shift;
  return defined($value) && $value =~ /^(?:1|true|yes)$/i;
}

1;
