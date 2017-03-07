def print_help = {
    log.info ''
    log.info 'GenomeArchitecture chip-seq pipeline'
    log.info '------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    chipseq.nf [dataset options] [genome options] [Other options]...'
    log.info ''
    log.info 'Dataset options:'
    log.info ''
    log.info '  Specify either --metadata or --datadir. Input/mock signals must have "input" as protein'
    log.info '  name. All input signals will be merged together and used as reference for all proteins.'
    log.info '  The supplied input signals must be compatible with each other and the experiments.'
    log.info ''
    log.info '    download datasets:'
    log.info '    --metadata METADATA_FILE     Comma separated file containing dataset information:'
    log.info '                                 Protein,Experiment,Replicate,read1URL[,read2URL]'
    log.info '    or:'
    log.info '    --datadir FASTQ_DIR          Filename format: Protein_Experiment_Replicate[_1,_2].fastq[.gz].'
    log.info ''
    log.info 'Reference genome options:'
    log.info ''
    log.info '  To use a previously-built BWA index, use --index.'
    log.info '  To automatically download and build, use --assembly with a supported genome (see below).'
    log.info '  To build the index from a local file, use --genome.'
    log.info '  To define the index output, use --index with --assembly/--genome. (Default: OUTPUT_DIR/genome).'
    log.info ''
    log.info '    provide bwa index path:'
    log.info '    --index INDEX_PATH           Path to a BWA index file. If not found, the index will be'
    log.info '                                 built from --assembly or --genome (if specified) and stored'
    log.info '                                 in INDEX_PATH.'
    log.info '    or automatic build (slow):'
    log.info '    --assembly GENOME_NAM        Supported genomes are hg19/hg38/mm9/mm10/dm3/dm6.'
    log.info '    or local file (slow):'
    log.info '    --genome GENOME_FASTA        Genome file in fasta format.'
    log.info ''
    log.info 'Other options:'
    log.info ''
    log.info '    --out OUTPUT_DIR             Output dir (Default: .).'
    log.info '    --confidence VALUE           Minimum ChIP enrichment confidence (Default: 0.99).'
    log.info '    --cpu N_CPUS                 Max CPU for mapping process (Default: 12).'
    log.info '    --colorspace                 Perform color space to nucleotide conversion (Default: false).'
    log.info '    --help                       Show this message and exit.'
    log.info ''
}


genome_urls = [ "hg19": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz",
                "hg38": "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz", 
                "mm9" : "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz",
                "mm10": "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz",
                "dm3" : "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz",
                "dm6" : "http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/chromFa.tar.gz"
              ]

// Set defaults.
params.out         = "./"
params.confidence  = 0.99
params.cpu         = 12
params.colorspace  = false
params.help        = false
params.metadata    = null
params.datadir     = null
params.index       = null
params.assembly    = null
params.genome      = null


if (params.help) {
   print_help()
   exit 0
}

if ((params.metadata && params.datadir) || !(params.metadata || params.datadir)) {
   log.info 'error: must specify either --metadata or --datadir.'
   log.info ''
   log.info 'use option --help for help.'
   exit 1
}

if (!(params.index || params.genome || params.assembly)) {
   log.info 'error: reference genome not provided.'
   log.info ''
   log.info 'use option --help for help.'
   exit 1
}

if (params.genome && params.assembly) {
   log.info 'error: specify either --assembly or --genome, not both.'
   log.info ''
   log.info 'use option --help for help.'
   exit 1
}

// 1. Parse metadata.

if (params.metadata) {
   Channel
   .fromPath(params.metadata)
   .splitCsv(strip: true)
   .into{dataset_list}

   process downloadDataset {
      // Process options
      tag "${info[0]}_${info[1]}_${info[2]}"
      publishDir path:"${params.out}/datasets/", mode:'symlink'
      // Cluster options
      cpus 1
      memory '4GB'

      input:
        val info from dataset_list
      output:
        set info, "*.fastq.gz" into fastq_files
      script:
        if (info.size() < 3) {
           log.info("format error in metadata file, ignoring line: ${info.join(',')}")
        }
        """
        fastq-dump --split-files --gzip -A ${info[2]}
        rm -f ~/ncbi/public/sra/${info[2]}.sra
        """
      /*
      if (info.size() == 4) {
        """
        wget ${info[3]} -O ${info[0]}_${info[1]}_${info[2]}.fastq.gz
        """
      } else if (info.size() == 5) {
        """
        wget ${info[3]} -O ${info[0]}_${info[1]}_${info[2]}_1.fastq.gz
        wget ${info[4]} -O ${info[0]}_${info[1]}_${info[2]}_2.fastq.gz
        """
      } else {
         log.info("format error in metadata file, ignoring line: ${info.join(',')}")
      }
      */
   }
} else {
   log.info 'sorry, --datadir option is not yet supported, please use --metadata.'
   log.info ''
   log.info 'use option --help for help.'
   exit 1
}


// 2. Find/Build BWA Index.
index_dir   = params.out + "/genome/"
index_name  = null
downl_fasta = Channel.create()
local_fasta = Channel.create()
bwa_index = Channel.create()
bwtnotfound = false

if (params.index) {
   // Check bwa index files.
   fasta_ref = file("${params.index}.bwt")
   index_ref = file(params.index)
   if (fasta_ref.isFile()) {
      bwa_index << params.index
   } else if (index_ref.getExtension() in ['bwt','amb','ann','pac','sa'] && index_ref.isFile()) {
      bwa_index << (index_ref.getParent().equals(null) ? './' : index_ref.getParent()) + index_ref.getBaseName()
   } else {
      bwtnotfound = true
   }
} 

if (!params.index || bwtnotfound) {
   if (params.assembly) {
      // Find assembly URL and trigger download.
      if (params.assembly in genome_urls.keySet()) {
         downl_fasta << genome_urls[params.assembly]
         index_name = params.assembly + ".fasta"
      } else {
         log.info 'error: specified assembly is not supported, use --genome instead.'
         log.info ''
         log.info 'use option --help for help.'
         exit 1
      }
   } else if (params.genome) {
      // Check file and trigger index build.
      if (file{params.genome}.isEmpty()) {
         local_fasta << params.genome
         index_name = file(params.genome).getName()
      } else {
         log.info 'error: genome reference file not found (${params.genome}).'
         log.info ''
         log.info 'use option --help for help.'
         exit 1
      }
   }
   if (params.index) {
      index_dir  = file(params.index).getParent()
      index_name = file(params.index).getName()
   }
}

process downloadReferenceGenome {
   // Process options
   //tag "${url}"
   // Cluster options
   cpus 1
   memory '4GB'

   input:
     val url from downl_fasta
   output:
     file "${index_name}.fasta" into local_fasta
   script:
     """
     wget ${url} -O ${index_name}.tar.gz
     tar -xvf ${index_name}.tar.gz
     cat *.fa > ${index_name}.fasta
     rm *.fa *.tar.gz
     """
}


process buildBWAindex {
   // Process options
   //tag "${index_path}"
   publishDir path: "${index_dir}", mode:'move'
   // Cluster options
   cpus 1
   memory '64GB'

   input:
     file "${index_name}" from local_fasta
   output:
     val index_path into bwa_index
     file "*.{bwt,amb,ann,pac,sa}" into index_files
   script:
     index_path = "${index_dir}/${index_name}"
     """
     bwa index ${index_name}
     """
}

// Map reads with BWA
process mapReads {
   // Process options
   tag "${info[2]}.bam"
   publishDir path:"${params.out}/mapped/", mode:'symlink'
   // Cluster options
   memory '32GB'
   clusterOptions "-pe smp ${params.cpu}"

  input:
    set info, files from fastq_files
    val index_path from bwa_index.first()
  output:
    set info, '*.bam' into bam_gene, bam_prot, bam_input
  script:
    filestr = files.size() == 2 ? "${files[0]} ${files[1]}" : "${files}"
    if (params.colorspace) {
       filestr = files.size() == 2 ? "<(colortobase ${files[0]}) <(${ctb_path} ${files[1]})" : "<(colortobase ${files})"
    }
    ref_file = (files.size() == 2 ? files[0] : files)
    if (ref_file.getExtension() == "gz") {
       seqlen_str = "zcat ${ref_file}"
    } else {
       seqlen_str = "cat ${ref_file}"
    }
    """
    OPTS_36='-k18 -B3 -O5 -T28'
    OPTS_26='-k17 -r1.3 -B2 -O4 -T22'
    OPTS_40=''
    OPTS_50=''
    SEQLEN=\$((\$(${seqlen_str} | head -2 | tail -1 | wc -c) - 1));
    if [ \$SEQLEN -le 30 ]; then \
      OPTS=\$OPTS_26; \
    elif [ \$SEQLEN -le 39 ]; then \
      OPTS=\$OPTS_36; \
    elif [ \$SEQLEN -le 46 ]; then \
      OPTS=\$OPTS_40; \
    else \
      OPTS=\$OPTS_50; \
    fi;
    echo \$SEQLEN
    echo \$OPTS
    bwa mem -t ${params.cpu} \$OPTS ${file(index_path)} ${filestr} | samtools view -bS - > ${info[2]}.bam
    """
}

// Filter files in groups (gene, experiment/lab, inputs)
bam_gene.filter{info,bam -> !(info[0] =~ /input/)}.map{info, bam -> tuple(info[0],bam)}.groupTuple().into{gene_signal}
bam_prot.filter{info,bam -> !(info[0] =~ /input/)}.map{info, bam -> tuple(info[0]+"_"+info[1],bam)}.groupTuple().into{prot_signal}
bam_input.filter{info,bam -> info[0] =~ /input/}.map{info, bam -> tuple(info[0],bam)}.groupTuple().into{input_signal}

// Make a combined signal array.
prot_signal.mix(gene_signal).into{chip_signal}

// Zerone files
process ZeroneDiscretization {
   // Process options
   tag "${file(chip).getName()}"
   publishDir path:"${params.out}/discretized/", mode:'move'
   // Cluster options
   cpus 1
   memory '16GB'

   input:
     set chip, sig_files from chip_signal
     set input, inp_files from input_signal.first()
   output:
     file "${chip}.01" into chip_out
     file "${chip}.bed" into chip_out_bed
   script:
     """
     zerone -c ${params.confidence} -0 ${inp_files.join(",")} -1 ${sig_files.join(",")} > ${chip}.01
     zerone -c ${params.confidence} -l -0 ${inp_files.join(",")} -1 ${sig_files.join(",")} > ${chip}.bed
     """
}
