# Current CEMT RNA-Seq RPKM computation

In: Genome wide coverage track (aligned to GRCh37lite)

Out: Normalized RNA-Seq track

* Selected transcripts model: ensembl69

Compute:

- Coverage over each exon (collapsed, no double counting of reads), excluding
  * Mitochondrial genome and genes coding ribosomal proteins
  * 0.5% of exons with highest coverage

Note RPKM[transcript] = [# of mapped reads]/([length of transcript]/10^3)/([total reads]/10^6)

N = [Total reads] = Total coverage over exonic regions not excluded

RPKM[exon] = Coverage[exon]/(length[exon]/10^3)/(N/10^6)

"RPKM[genomic_position]" = Coverage[genomic_position]/(1/10^3)/(N/10^6)

Define normalized RNA-Seq track as Coverage_Track/(10^3/N)

(Coverage tracks reporting one uniform RPKM value for the gene may also reported)



 

Key observation: Integrating the normalized track over an exon yields RPKM[exon]
