# IsoSeq-3.0 : Generating full-length cDNA sequences for transcriptomics

The isoform sequencing (IsoSeq) application generates full-length cDNA sequences — from the 5’ end of transcripts to the poly-A tail — eliminating the need for transcriptome reconstruction using isoform-inference algorithms. The Iso-Seq method generates accurate information about alternatively spliced exons and transcriptional start sites. It also delivers information about poly-adenylation sites for transcripts up to 10 kb in length across the full complement of isoforms within targeted genes or the entire transcriptome.

Table of contents
=================

  * [Overview](#overview)
  * [Manual](#manual)
    * [Running with SMRTLink](#running-with-smrtlink)
    * [Running on the Command-Line](#running-on-the-command-line)
    * [Running on the Command-Line with PBSMRTPipe](#running-on-the-command-line-with-pbsmrtpipe)
  * [Options](#options)
    * [Classify](#classify-options)
    * [Cluster](#cluster-options)
    * [Subset](#subset-options)
  * [Files](#files)
    * [Classify](#classify-files)
    * [Cluster](#cluster-files)
  * [Algorithms](#algorithms)
  * [Glossary](#glossary)



## Overview
Analyses are performed in three stages, CCS, Classify and Cluster. For analyses performed on the command-line, there is an optional tool, Subset, for subsetting the IsoSeq results.
* __CCS__
  * CCS is the first stage of an IsoSeq analysis. CCS builds circular consensus sequences from your subreads. More information about CCS is available here: https://github.com/PacificBiosciences/pbccs/blob/master/README.md
* __Classify__
  * Classify is the second stage of an IsoSeq analysis. The key output of Classify is a file of full-length non-chimeric reads, and a file of non-full length reads. The key input of Classify is the circular consensus sequences generated from CCS. Classify will identify and remove polyA/T tails, remove primers, and identify read strandedness. Classify also removes artificial concatemers, but does not remove PCR chimeras. 
* __Cluster__
  * Cluster is the third stage of an IsoSeq analysis. The key outputs of Cluster is a file of polished, high-quality consensus sequences, and a file of polished, low-quality consensus sequences. The key input of clustering is the file of full-length non-chimeric reads, and a file of non-full length reads outputted by Classify. 
* __Subset__
  * Subset is an optional program which can be used to subset the output files for particular classes of sequences, such as non-chimeric reads, or non-full-length reads.

##Manual

There are three ways to run IsoSeq: Using SMRTLink, on the command-line, and on the command-line using pbsmrtpipe so that you can run the whole IsoSeq job with one command given to pbsmrtpipe. 

##Running with SMRTLink

##Running on the Command-Line

Without PBSMRTPipe, the analysis is performed in 3 steps:

1. Run CCS on your subreads, generating a CCS BAM file. Then generate an XML from the BAM file.
2. Run Classify on your CCSs with the XML as input, generating a FASTA of annotated sequences.
3. Run Cluster on the FASTA produced by Classify, generating polished isoforms. 

__Step 1. CCS__

First convert your subreads to circular consensus sequences. You can do this with the command:

     ccs ccs.bam subreads.bam

Where ccs.bam is where the CCSs will be output, and subreads.bam is the file containing your subreads. 
Next, you will generate an XML file from your CCSs. You can do this with the commmand:

     dataset create --type ConsensusReadSet ccs.xml ccs.bam

Where ccs.xml is the name of the xml file you are generating and ccs.bam is the name of the bam file you generated previously using the ccs command. 

__Step 2. Classify__

Classify can be run at the command line as follows:

     pbtranscript classify [OPTIONS] ccs.xml classified.fasta

Where ccs.xml is the xml file you generated in Step 1, and classified.fasta is your output file. 

__Step 3. Cluster__

Cluster can be run at the command line as follows:

     pbtranscript cluster [OPTIONS] flnc_fa consensusFa

__Optional Step__

Once Cluster has run, you can further subset your results using Subset. Subset can be run with the command:

    pbtranscript subset [OPTIONS] isoseq_draft.fasta isoseq_subset.fasta

Where isoseq_draft.fasta is the input FASTA and isoseq_subset.fasta is the output FASTA. 

##Running on the Command-Line with PBSMRTPipe

For the user who would like a simplified experience, pbsmrtpipe offers a way to run isoseq with a single command. The caveat to this approach is that only a few options are available to the user. 
run isoseq with pbsmrtpipe, first load the smrtpipe module.

```
 module load smrtanalysis/3.0.1-current
```
Now create an xml file from your subreads.

```
dataset create --type SubreadSet subreads subreads.bam
```
This will create a file called subreads.xml. Now create a global options xml file and an isoseq options xml file.

```
 pbsmrtpipe show-workflow-options -o global_options.xml
 pbsmrtpipe show-template-details pbsmrtpipe.pipelines.sa3_ds_isoseq -o isoseq_options.xml
```

The options you may modify are now contained in the files global_options.xml and isoseq_options.xml. An entry in these files looks like this:

```
 <option id="pbtranscript.task_options.min_seq_len">
            <value>300</value>
        </option>
```

And you can modify them using your favorite text editor, such as vim.
Once you have set your options, you are ready to run isoseq via pbsmrtpipe:

```
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_isoseq -e eid_subread:subreads.xml --preset-xml=isoseq_options.xml --preset-xml=global_options.xml
```

## Options
## Classify Options

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
| Summary | -summary SUMMARY_FN.txt | TXT file to output classify summary (default: out.classify_summary.txt) |
| Primers File | -p PRIMERFN, --primer PRIMERFN  | Primer fasta file (default: primers.fa) |
| Primers Report | --report PRIMERREPORTFN  | CSV file of primer info. Contains the same info found in the description lines of the output FASTA (default: out.primer_info.csv) |
| CPUs | --cpus CPUS  | Number of CPUs to run HMMER (default: 8) |



|      Chimera-detection Arguments     |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Minimum Sequence Length  | --min_seq_len MIN_SEQ_LEN   | Minimum sequence length to output (default: 300) |
| Minimum PHMMER Score  | --min_score MIN_SCORE   | Minimum phmmer score for primer hit (default: 10) |
| Non-Full-Length Chimeras  | --detect_chimera_nfl   | Detect chimeric reads among non-full-length reads. Non-full-length non-chimeric/chimeric reads will saved to outDir/nflnc.fasta and outDir/nflc.fasta. |

|      Read-Extraction Arguments     |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| Ignore polyA  | --ignore_polyA   | FL does not require polyA tail (default: turned off) |


## Cluster Options

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
| Report  | --report REPORT_FN | NO DESCRIPTION (yli) |
| Pickle???  | --pickle_fn PICKLE_FN | NO DESCRIPTION (yli) |

|           ICE Arguments           |     Example      |  Explanation      |
| -------------------------- | --------------------------- | ----------------- |
| cDNA  | --cDNA_size {under1k,between1k2k,between2k3k,above3k} | Estimated cDNA size. |
| Quiver  | --quiver | Call quiver to polish consensus isoforms using non-full-length non-chimeric CCS reads. |
| Finer Quiver  | -h, --help | Use finer classes of QV information from CCS input instead of a single QV from FASTQ. This option is slower and consumes more memory. |

|           SGE environment Arguments          |     Example      |  Explanation      |
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


## Subset Options

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

## Files
## Classify Files
__Output FASTA (out.fasta)__
Reads from Classify look like this:

```
>m140121_100730_42141_c100626750070000001823119808061462_s1_p0/119/30_1067_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=1067;threeend=1096;primer=1;chimera=0
ATAAGACGACGCTATATG
 ```
These lines have the format: 
```
<movie_name>/<ZMW>/<start>_<end>_CCS INFO
```
The info fields are:
* strand: either + or -,  whether a read is forward or reverse-complement cdna,
* fiveseen: whether or not 5' prime is seen in this read, 1 yes, 0 no
* polyAseen: whether or not poly A tail is seen, 1 yes, 0 no
* threeseen: whether or not 3' prime is seen, 1 yes, 0 no
* fiveend: start position of 5'
* threeend: start position of 3' in read
* polyAend: start position of polyA in read
* primer: index of primer seen in this read (remember  primer fasta file >F0 xxxxx >R0 xxxxx >F1 xxxxx >R1 xxxx)
* chimera: whether or not this read is classified as  a chimeric cdna

__Summary (out.classify_summary.txt)__ 
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

##Cluster Files

__Output Isoforms__
todo (yli)

__Summary__
todo (yli)

__Report__
todo (yli)



## Algorithms

__CCS__

See https://github.com/PacificBiosciences/pbccs/blob/master/README.md

__Classify__

todo (yli)

__Cluster__

todo (yli)

__HMMER__

HMMER is used for searching sequence databases for homologs of protein sequences, and for making protein sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs).

__ICE__

todo (yli)

__Quiver__

maybe explain briefly how Quiver is used in isoseq, a more general explanation of Quiver can exist elsewhere (yli)

## Glossary
* __Chimera__
  * todo (yli)
* __Concatemer__
  * todo (yli)
