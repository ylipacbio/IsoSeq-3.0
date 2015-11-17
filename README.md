# IsoSeq-3.0 : Generating full-length cDNA sequences for transcriptomics

The isoform sequencing (IsoSeq) application generates full-length cDNA sequences — from the 5’ end of transcripts to the poly-A tail — eliminating the need for transcriptome reconstruction using isoform-inference algorithms. The Iso-Seq method generates accurate information about alternatively spliced exons and transcriptional start sites. It also delivers information about poly-adenylation sites for transcripts up to 10 kb in length across the full complement of isoforms within targeted genes or the entire transcriptome.

## Command-Line Overview
## Classify, Cluster, Subset
Analyses are performed using three tools
* __Classify__
  * Classify is the first program to be run when performing an IsoSeq analysis. The key output of Classify is a file of full-length non-chimeric reads, and a file of non-full length reads. The key input of Classify is the set of subreads produced from running CCS on the subreads from your PacBio instrument. Classify will identify and remove polyA/T tails, remove primers, and identify read strandedness. Classify also removes artificial concatemers, but does not remove PCR chimeras. 
* __Cluster__
  * Cluster is the second program to be run when performing an IsoSeq analysis. The key outputs of Cluster is a file of polished, high-quality consensus sequences, and a file of polished, low-quality consensus sequences. The key input of clustering is the file of full-length non-chimeric reads, and a file of non-full length reads outputted by Classify. 
* __Subset__
  * Subset is an optional program which can be used to subset the output files for particular classes of sequences, such as non-chimeric reads, or non-full-length reads. 

## Files

__Classify Output FASTA__
Reads from classify look like this

```
     m140121_100730_42141_c100626750070000001823119808061462_s1_p0/119/30_1067_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=1067;threeend=1096;primer=1;chimera=0
     
 ```
The first field has the format: 
```
<movie_name>/<ZMW>/<start>_<end>_CCS INFO
```
The info fields are:
* strand: either + or -,  whether a read is forward or reverse-complement cdna,
* fiveseen: whether or not 5' prime is seen in this read, 1 yes, 0 no
* polyAseen: whether or not poly A tail is seen
* threeseen: whether or not 3' prime is seen
* fiveend: start position of 5'
* threeend: start position of 3' in read
* polyAend: start position of polyA in read
* primer: index of primer seen in this read (remember  primer fasta file >F0 xxxxx >R0 xxxxx >F1 xxxxx >R1 xxxx)
* chimera: whether or not this read is classified as  a chimeric cdna



__Classify Summary File__ 
This file contains the following statistics:
* Number of reads of insert
* Number of five prime reads
* Number of three prime reads
* Number of poly-A reads
* Number of filtered short reads
* Number of non-full-length reads
* Number of full-length reads
* Number of full-length non-chimeric reads
* Average full-length non-chimeric read length




## Command-Line Manual
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
| Summary | -summary SUMMARY_FN.txt | TXT file to output classify summary (default: classify_summary.txt) |
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

The input file and the output file must be fasta format. An example command would be:

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
| CCS QVs FOFN  | --ccs_fofn CCS_FOFN | A FOFN of ccs.h5 or ccs.bam (e.g., ccs.fofn), which contain quality values of consensus (CCS) reads. If not given, assume there is no QV information available. |
| Reads QVs FOFN |  --bas_fofn BAS_FOFN  | A FOFN of bax/bas.h5 or bam files (e.g., input.fofn), which contain quality values of raw reads and subreads |
| Output Directory  | -d ROOT_DIR, --outDir ROOT_DIR | Directory to store temporary and output cluster files.(default: output/) |
| Temp Directory  | --tmp_dir TMP_DIR | Directory to store temporary files.(default, write to root_dir/tmp.). |
| Summary  | --summary SUMMARY_FN | TXT file to output cluster summary (default: my.cluster_summary.txt) |
| Report  | --report REPORT_FN | NO DESCRIPTION |
| Pickle???  | --pickle_fn PICKLE_FN | NO DESCRIPTION |

|           ICE Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| cDNA  | --cDNA_size {under1k,between1k2k,between2k3k,above3k} | Estimated cDNA size. |
| Quiver  | --quiver | Call quiver to polish consensus isoforms using non-full-length non-chimeric CCS reads. |
| Finer Quiver  | -h, --help | Use finer classes of QV information from CCS input instead of a single QV from FASTQ. This option is slower and consumes more memory. |

|           ICE Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Run SGE  | --use_sge | Instructs Cluster to use SGE |
| Maximum SGE Jobs  | --max_sge_jobs MAX_SGE_JOBS | The maximum number of jobs that will be submitted to SGE concurrently. |
| SGE Job ID  | --unique_id UNIQUE_ID | Unique ID for submitting SGE jobs. |
| BLASR Cores  | --blasr_nproc BLASR_NPROC | Number of cores for each BLASR job. |
| Quiver CPUs  | --quiver_nproc QUIVER_NPROC | Number of CPUs each quiver job uses. |

|           IceQuiver High QV/Low QV Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Minimum Quiver Accuracy  | --hq_quiver_min_accuracy HQ_QUIVER_MIN_ACCURACY | Minimum allowed quiver accuracy to classify an isoform as hiqh-quality. |
| Trim QVs 5'  | --qv_trim_5 QV_TRIM_5 | Ignore QV of n bases in the 5' end. |
| Trim QVs 3'  | --qv_trim_3 QV_TRIM_3 | Ignore QV of n bases in the 3' end. |
| High-Quality Isoforms FASTA  | --hq_isoforms_fa HQ_ISOFORMS_FA | Quiver polished, high quality isoforms in fasta, default: root_dir/output/all_quivered_hq.fa |
| High-Quality Isoforms FASTQ  | --hq_isoforms_fq HQ_ISOFORMS_FQ | Quiver polished, high quality isoforms in fastq, default: root_dir/output/all_quivered_hq.fq |
| Low-Quality Isoforms FASTA  | --lq_isoforms_fa LQ_ISOFORMS_FA | Quiver polished, low quality isoforms in fasta, default: root_dir/output/all_quivered_lq.fa |
| Low-Quality Isoforms FASTQ  | --lq_isoforms_fq LQ_ISOFORMS_FQ | Quiver polished, low quality isoforms in fastq, default: root_dir/output/all_quivered_lq.fq |


## Subset

Subset can be run at the command line as follows:

     pbtranscript subset [OPTIONS] readsFN outFN

The input file and the output file must be fasta format. An example command would be:

    pbtranscript subset [OPTIONS] isoseq_draft.fasta isoseq_subset.fasta
    
Cluster can be run with a variety of options described below.  

|           Positional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Input Sequences  | isoseq_draft.fasta  | Input fasta file (usually isoseq_draft.fasta) |
| Output Sequences | isoseq_subset.fasta | Output fasta file |

|           Optional Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Help  | -h, --help | This prints the help message |
| Output Full-length  | --FL | Reads to output must be Full-Length, with 3' primer and 5' primer and polyA tail seen. |
| Output Non-Full-length  | --nonFL | Reads to output must be Non-Full-Length reads. |
| Output Non-Chimeric  | --nonChimeric | Reads to output must be non-chimeric reads. |
| Output Read-Length  | --printReadLengthOnly | Only print read lengths, no read names and sequences. |
| Ignore polyA Tails  |  --ignore_polyA | FL does not require polyA tail (default: turned off) |

