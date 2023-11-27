=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 SpliceVault

=head1 SYNOPSIS

 mv SpliceVault.pm ~/.vep/Plugins

 ./vep -i variations.vcf --plugin SpliceVault,file=/path/to/SpliceVault_data.tsv.gz

 # Stringely select predicted loss-of-function (pLoF) splicing variants
 ./filter_vep -i variant_effect_output.txt --filter "SPLICEVAULT_OUT_OF_FRAME_EVENTS >= 3"

=head1 DESCRIPTION

 A VEP plugin that retrieves SpliceVault data to predict exon-skipping events
 and activated cryptic splice sites based on the most common mis-splicing events
 around a splice site.

 This plugin returns the most common (top 4) variant-associated mis-splicing
 events for variants that overlap a near-splice-site region (-17:+3 for splice
 acceptor sites and -3:+8 for splice donor sites) based on SpliceVault data.
 Each event includes the following information:
   - Type: exon skipping (ES), cryptic donor (CD) or cryptic acceptor (CA)
   - Description: exons skiped (ES events) or cryptic distance to the annotated
     splice site (CD/CA events)
   - Percent of supporting samples: samples supporting the event over total
     samples where splicing occurs in that site (in percentage)
   - Frameshift: inframe or out-of-frame event

 Please cite the SpliceVault publication alongside the VEP if you use this
 resource: https://pubmed.ncbi.nlm.nih.gov/36747048/

 The tabix utility must be installed in your path to use this plugin. The
 SpliceVault TSV and respective index (TBI) for GRCh38 can be downloaded from
 https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data.tsv.gz

 To filter results, please use filter_vep with the output file or standard
 output. Documentation on filter_vep is available at:
 https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

=cut

package SpliceVault;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();
  
  # Check files in arguments
  my @files; 
  my $params = $self->params_to_hash();

  for my $key (keys %{$params}) {
    push @files, $params->{$key};
  }

  die "ERROR: add Tabix-indexed file, for example:\n" .
    "  --plugin SpliceVault,file=/path/to/SpliceVault_data.tsv.gz \n"
    unless @files > 0;
  $self->add_file($_) for @files;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    SpliceVault_top_events          => 'List of SpliceVault top events with the following colon-separated fields per event: rank:type:transcript_impact:percent_of_supporting_samples:frame',
    SpliceVault_site_sample_count   => 'Number of SpliceVault annotated splicing samples',
    SpliceVault_site_max_depth      => 'Maximum number of reads seen in any one sample representing annotated splicing in GTEx',
    SpliceVault_out_of_frame_events => 'Fraction of SpliceVault top events that cause a frameshift; as per SpliceVault article, sites with 3/4 or more in-frame events are likely to be splice-rescued and not LoF',
    SpliceVault_site_pos            => 'Genomic position of splice-site predicted to be lost by SpliceAI',
    SpliceVault_site_type           => 'Type of splice-site predicted to be lost by SpliceAI',
  }
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $allele = $tva->base_variation_feature->alt_alleles;

  my @data = @{ $self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end}) };

  foreach (@data) {
    next unless $tva->transcript->stable_id eq $_->{feature};

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $vf->ref_allele_string,
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    next unless (@$matches);

    if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
      my $top_events = $_->{result}->{SpliceVault_top_events};
      my $res = {};
      for my $event (@$top_events) {
        my ($rank, $type, $impact, $samples, $frame) = split(/:/, $event);
        $res->{$rank}->{'site_type'}                     = $type;
        $res->{$rank}->{'transcript_impact'}             = $impact;
        $res->{$rank}->{'percent_of_supporting_samples'} = $samples;
        $res->{$rank}->{'frame'}                         = $frame;
      }
      $_->{result}->{SpliceVault_top_events} = $res;
      return { 'SpliceVault' => $_->{result} };
    } else {
      return $_->{result};
    }
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $start, $alt, $feature, $type, $site, $out_of_frame, $top_events, $count, $max_depth) = split /\t/, $line;

  $top_events =~ s/ /_/g;
  $top_events =~ s/;/:/g;

  my $res = {
    feature => $feature,
    alt     => $alt,
    start   => $start,
    result  => {
      SpliceVault_site_sample_count   => $count,
      SpliceVault_site_max_depth      => $max_depth,
      SpliceVault_out_of_frame_events => $out_of_frame,
      SpliceVault_top_events          => [ split(/\|/, $top_events) ],
      SpliceVault_site_pos            => $site,
      SpliceVault_site_type           => $type,
    }
  };
  return $res;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
