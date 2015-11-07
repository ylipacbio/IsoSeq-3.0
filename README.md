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
| Output Non-Full-Length  | --nfl NFL_FA.fasta | Outputs non-full-length reads in fasta |

|           HMMER Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| HMMER Directory | -d OUTDIR, --outDir OUTDIR  | Directory to store HMMER output (default: output/) |
| Summary | -summary SUMMARY_FN.txt | TXT file to output classsify summary (default: classify_summary.txt) |
| Primers File | -p PRIMERFN, --primer PRIMERFN  | Primer fasta file (default: primers.fa) |
| Primers Report | --report PRIMERREPORTFN  | CSV file to output primer info (default: .primer_info.csv) |
| CPUs | --cpus CPUS  | Number of CPUs to run HMMER (default: 8) |



|      Chimera-detection Arguments     |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Minimum Sequence Length  | --min_seq_len MIN_SEQ_LEN   | Minimum sequence length to output (default: 300) |
| Minimum PHMMER Score  | --min_score MIN_SCORE   | Minimum phmmer score for primer hit (default: 10) |
| Non-Full-Length Chimeras  | --detect_chimera_nfl   | Detect chimeric reads among non-full-length reads. Non-full-length non-chimeric/chimeric reads will saved to outDir/nflnc.fasta and outDir/nflc.fasta. |

|      Read-Extraction Arguments     |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Ignore polyA  | --ignore_polyA   | FL does not require polyA tail (default: turned off) |


## Cluster

Cluster can be run at the command line as follows:

     pbtranscript cluster [OPTIONS] flnc_fa consensusFa

The input file can be fasta or bam format, and the output file must be fasta format. An example command would be:

    pbtranscript cluster [OPTIONS] isoseq_flnc.fasta output.fasta
    
Cluster can be run with a variety of options described below.  

|           Positional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Input Reads  | isoseq_flnc.fasta  | Input full-length non-chimeric reads in fasta format, used for clustering consensus isoforms |
| Output Isoforms | out.fasta | Output predicted (unpolished) consensus isoforms in fasta file. |

|           Optional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Help  | -h, --help | This prints the help message |
| Input Non-Full-Length  | --nfl_fa NFL_FA.fasta | Input non-full-length reads in fasta format, used for polishing consensus isoforms, e.g., isoseq_nfl.fasta |
| Quality Values FOFN  | --ccs_fofn CCS_FOFN | A FOFN of ccs.h5 or ccs.bam (e.g., ccs.fofn), which contain quality values of consensus (CCS) reads. If not given, assume there is no QV information available. |
