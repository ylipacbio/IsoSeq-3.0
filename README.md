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

|           Positional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| readsFN  | ccs.bam  | This is the second-to-last argument and it names the file containing the input sequence data in BAM or FASTA format |
| outReadsFN | out.fasta | This is the last argument and it names the file containing the input sequence data in BAM or FASTA format |

|           Optional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Help  | -h, --help | This prints the help message |
| Full-Length Non-Chimeric  | --flnc FLNC_FA.fasta | Outputs full-length non-chimeric reads in fasta |
| Non-Full-Length  | --nfl NFL_FA.fasta | Outputs non-full-length reads in fasta |

|           HMMER Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| HMMER Directory | -d OUTDIR, --outDir OUTDIR  | Directory to store HMMER output (default: output/) |
| Summary | -summary SUMMARY_FN.txt | TXT file to output classsify summary (default: classify_summary.txt) |
| Primers File | --report PRIMERREPORTFN.csv  | CSV file to output primer info (default: .primer_info.csv) |
| Primers Report | --cpus CPUS  | Number of CPUs to run HMMER (default: 8) |
| CPUs | --cpus CPUS  | Number of CPUs to run HMMER (default: 8) |



|           Chimera-detection Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Output File Name           | myResult.bam                | This argument is the first argument that comes after the named arguments shown below. |
