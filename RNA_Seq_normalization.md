# Current CEMT RNA-Seq RPKM computation

In: Genome wide coverage track (aligned to transcriptome and repositioned to GRCh37lite)

Out: Normalized RNA-Seq track

* Selected transcripts model: ensembl69

Compute:

- Read coverage over each exon (exons collapsed, no double counting of reads) normalized by readlength, excluding
  * Mitochondrial genome and genes coding ribosomal proteins
  * 0.5% of exons with highest coverage

Note RPKM[transcript] = [# of mapped reads]/([length of transcript]/10^3)/([total reads]/10^6)

N = [total reads] = Total normalized coverage over exonic regions not excluded by above criterion

RPKM[exon] = normalized_coverage[exon]/(length[exon]/10^3)/(N/10^6)

"RPKM[genomic_position]" = normalized_coverage[genomic_position]/(1/10^3)/(N/10^6)

Define normalized RNA-Seq track as normalized_coverage_track/(10^9/N)

(Coverage tracks reporting one uniform RPKM value for the gene may also reported)



 

Key observation: Integrating the normalized track over an exon yields RPKM[exon]
