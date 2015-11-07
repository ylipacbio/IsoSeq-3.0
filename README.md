# IsoSeq-3.0 : Generating full-length cDNA sequences for transcriptomics

The isoform sequencing (Iso-Seq) application generates full-length cDNA sequences — from the 5’ end of transcripts to the poly-A tail — eliminating the need for transcriptome reconstruction using isoform-inference algorithms. The Iso-Seq method generates accurate information about alternatively spliced exons and transcriptional start sites. It also delivers information about poly-adenylation sites for transcripts up to 10 kb in length across the full complement of isoforms within targeted genes or the entire transcriptome.

## Command-Line Overview
Analyses are performed using three tools
* Classify
  * description of Classify
* Cluster
  * description of Cluster
* Subset
  * description of Subset

## Classify

Classify can be run at the command line as follows:

     pbtranscript classify [OPTIONS] input_file output_file

The input file can be fasta or bam format, and the output file must be fasta format. An example command would be:

    pbtranscript classify [OPTIONS] ccs.bam output.fasta
    
Classify can be run with a variety of options described below.  

|           Option           |     Example (Defaults)      |                                                                                                                                                                                                                                                                                                    Explanation                                                                                                                                                                                                                                                                                                     |
| -------------------------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Output File Name           | myResult.bam                | This argument is the first argument that comes after the named arguments shown below.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
