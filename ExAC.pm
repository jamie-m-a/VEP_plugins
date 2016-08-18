=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 ExAC

=head1 SYNOPSIS

 mv ExAC.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,,AN
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC,AN



=head1 DESCRIPTION

 A VEP plugin that retrieves ExAC allele frequencies.
 
 Visit ftp://ftp.broadinstitute.org/pub/ExAC_release/current to download the latest ExAC VCF.
 
 Note that the currently available version of the ExAC data file (0.3) is only available
 on the GRCh37 assembly; therefore it can only be used with this plugin when using the
 VEP on GRCh37. See http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly
 
 The tabix utility must be installed in your path to use this plugin.

 The plugin takes 3 command line arguments. Second and third arguments are not mandatory. If AC specified as second
 argument Allele counts per population will be included in output. If AN specified as third argument Allele specific
 chromosome counts will be included in output.


=cut

package ExAC;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  
  # get ExAC file
  my $file = $self->params->[0];

  # get AC,AN options
  if (exists($self->params->[1]) && $self->params->[1] eq 'AC'){
    $self->{display_ac} = 1;
  }
  else {
    $self->{display_ac} = 0;
  }

  if (exists($self->params->[2]) && $self->params->[2] eq 'AN'){
    $self->{display_an} = 1;
  }
  else {
    $self->{display_an} = 0;
  }

  # remote files?
  if($file =~ /tp\:\/\//) {
    my $remote_test = `tabix -f $file 1:1-1 2>&1`;
    if($remote_test && $remote_test !~ /get_local_version/) {
      die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
    }
  }

  # check files exist
  else {
    die "ERROR: ExAC file $file not found; you can download it from ftp://ftp.broadinstitute.org/pub/ExAC_release/current\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }
  
  $self->{file} = $file;
  
  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  
  if(!exists($self->{header_info})) {
    open IN, "tabix -f -h ".$self->{file}." 1:1-1 |";
    
    my %headers = ();
    my @lines = <IN>;
    
    while(my $line = shift @lines) {
      if($line =~ /ID\=AC(\_[A-Zdj]+)?\,.*\"(.+)\"/) {
        my ($pop, $desc) = ($1, $2);
        
        $desc =~ s/Counts?/frequency/i;
        $pop ||= '';
        
        my $field_name = 'ExAC_AF'.$pop;
        $headers{$field_name} = 'ExAC '.$desc;

        if ($self->{display_ac}){
          $field_name = 'ExAC_AC'.$pop;
          $headers{$field_name} = 'ExAC'.$pop.' Allele count';
        }
        if ($self->{display_an}){
          $field_name = 'ExAC_AN'.$pop;
          $headers{$field_name} = 'ExAC'.$pop.' Allele number';
        }

        # store this header on self
        push @{$self->{headers}}, 'AC'.$pop;
      }
    }
    
    close IN;
    
    die "ERROR: No valid headers found in ExAC VCF file\n" unless scalar keys %headers;
    
    $self->{header_info} = \%headers;
  }
  
  return $self->{header_info};
}

sub run {
  my ($self, $tva) = @_;

  # make sure headers have been loaded
  $self->get_header_info();

  my $vf = $tva->variation_feature;
  
  # get allele, reverse comp if needed
  my $allele;
  
  $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  # adjust coords to account for VCF-like storage of indels
  my ($s, $e) = ($vf->{start} - 1, $vf->{end} + 1);
  
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $s, $e);
  
  # clear cache if it looks like the coords are the same
  # but allele type is different
  delete $self->{cache} if
    defined($self->{cache}->{$pos_string}) &&
    scalar keys %{$self->{cache}->{$pos_string}} &&
    !defined($self->{cache}->{$pos_string}->{$allele});
  
  my $data = {};
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    $data = $self->{cache}->{$pos_string};
  }
  
  # read from file
  else {
    open TABIX, sprintf("tabix -f %s %s |", $self->{file}, $pos_string);
    
    while(<TABIX>) {
      chomp;
      s/\r$//g;
      
      # parse VCF line into a VariationFeature object
      my ($vcf_vf) = @{parse_line({format => 'vcf'}, $_)};
      
      # check parsed OK
      next unless $vcf_vf && $vcf_vf->isa('Bio::EnsEMBL::Variation::VariationFeature');
      
      # compare coords
      next unless $vcf_vf->{start} == $vf->{start} && $vcf_vf->{end} == $vf->{end};
      
      # get alleles, shift off reference
      my @vcf_alleles = split /\//, $vcf_vf->allele_string;
      my $ref_allele = shift @vcf_alleles;
      
      # iterate over required headers
      HEADER:
      foreach my $h(@{$self->{headers} || []}) {
        my $total_ac = 0;
        
        if(/$h\=([0-9\,]+)/) {
          
          # grab AC
          my @ac = split /\,/, $1;
          next unless scalar @ac == scalar @vcf_alleles;
          
          # now sed header to get AN
          my $anh = $h;
          $anh =~ s/AC/AN/;
          
          my $afh = $h;
          $afh =~ s/AC/AF/;

          # get AC from header
          my $ach = $h;

          if(/$anh\=([0-9\,]+)/) {
            
            # grab AN
            my @an = split /\,/, $1;
            next unless @an;
            my $an;

            foreach my $a(@vcf_alleles) {
              my $ac = shift @ac;
              $an = shift @an if @an;
              next HEADER unless $an;

              $total_ac += $ac;
              if ($self->{display_ac}){
                $data->{$a}->{'ExAC_'.$ach} = $ac;
              }
              if ($self->{display_an}){
                $data->{$a}->{'ExAC_'.$anh} = $an;
              }

              $data->{$a}->{'ExAC_'.$afh} = sprintf("%.3g", $ac / $an);
            }
            
            # use total to get ref allele freq
            if ($self->{display_ac}){
             $data->{$ref_allele}->{'ExAC_'.$ach} = $total_ac;
            }
            if ($self->{display_an}){
              $data->{$ref_allele}->{'ExAC_'.$anh} = $an;
            }
            $data->{$ref_allele}->{'ExAC_'.$afh} = sprintf("%.3g", 1 - ($total_ac / $an));
          }
        }
      }
    }
    
    close TABIX;
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => $data};
  
  return defined($data->{$allele}) ? $data->{$allele} : {};
}

1;

