#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_block = params.monochrome_logs ? '' : "\033[3m";
  c_ul = params.monochrome_logs ? '' : "\033[4m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_bul = c_bold + c_ul;
  is_viruses = (params.taxids == 10239) ? " (Viruses)" : ""
  log.info"""
  =${c_dim}=================================================================${c_reset}
  ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}==================================================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bul}Usage:${c_reset}
  The typical command for running the pipeline is as follows:
  
    nextflow run ${workflow.manifest.name} \\
      ${c_red}--reads "${params.reads}"${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      -profile singularity # Recommended to run workflow with Singularity

  The above ${c_bul}assumes${c_reset} that you have a ${c_cyan}Centrifuge DB${c_reset} and ${c_purple}Kraken2 DB${c_reset} located at
  ${c_cyan}/opt/DB/centrifuge/nt-2018-03-03/nt${c_reset} and ${c_purple}/opt/DB/kraken2/standard2${c_reset}, 
  respectively, ${c_bul}OR${c_reset} that you have set ${c_cyan}\$CENTRIFUGE_DB${c_reset} and ${c_purple}\$KRAKEN2_DB${c_reset} env 
  variables. It also assumes that you have ${c_yellow+c_bul}Singularity${c_reset} installed on your
  local machine and will automatically pull and use the Singularity image for
  this workflow from Singularity-Hub.org.

  ${c_yellow+c_bold+c_block}NOTE:${c_yellow} For best results, please ensure you have ${c_bul}Singularity${c_yellow} installed prior to running this workflow.${c_dim}(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)${c_reset}


  ${c_bul}Mandatory Options:${c_reset}
    ${c_red}--reads${c_reset}   Input reads directory and pattern (default: ${c_red}"${params.reads}"${c_reset})

  ${c_bul}Read Quality Filtering Options:${c_reset}
    --fastp_min_base_quality            Minimum base quality score for fastp read
                                        filtering. (default: ${params.fastp_min_base_quality})
    --fastp_max_percent_low_qual_base   Filter reads with more than X% low 
                                        quality bases as defined by 
                                        `--fastp_min_base_quality` (0-100) 
                                        (default: ${params.fastp_max_percent_low_qual_base})

  ${c_bul}Taxonomic Classification Options:${c_reset}
    ${c_cyan}--centrifuge_db${c_reset}   Path to Centrifuge DB and prefix. If not specified, will 
                      try to get from \$CENTRIFUGE_DB env variable or see if
                      "/opt/DB/centrifuge/nt-2018-03-03/nt" exists.
                      (default: ${c_cyan}${params.centrifuge_db}${c_reset})
    ${c_purple}--kraken2_db${c_reset}      Path to Kraken2 DB directory. . If not specified, will 
                      try to get from \$KRAKEN2_DB env variable or see if
                      "/opt/DB/kraken2/standard2" exists.
                      (default: ${c_purple}${params.kraken2_db}${c_reset})
    --taxids          Taxonomic IDs to filter reads by. Multiple taxids should
                      be delimited by commas (`--taxids 1,2,3`). To disable 
                      filtering of reads based on taxids, do not provide a
                      value for the `--taxids` argument:
                      `nextflow run ... --taxids --reads ...`
                      (default: ${params.taxids}${is_viruses})
    --exclude_unclassified_reads  Exclude unclassified reads from taxonomic
                                  classification filtered reads (default: false)

  ${c_bul}BLAST Options:${c_reset}
    --blastn_db        Nucleotide BLAST database (e.g. /opt/DB/blast/nt if NCBI nt DB downloaded to /opt/DB/blast/)
    --blastn_taxids    File containing taxids to restrict nucleotide BLAST search (default: ${params.blastn_taxids}).
                       If you do not want to restrict the BLAST search, set this parameter to ${c_yellow}--blastn_taxids ''${c_reset}
    --blastn_qcov_hsp_perc  BLASTN query HSP percent coverage threshold (default: ${params.blastn_qcov_hsp_perc})

  ${c_bul}De Novo Assembly Options:${c_reset}
    --unicycler_mode  Unicycler assembly mode (default: ${params.unicycler_mode})
    --shovill_trim    Trim reads with Trimmomatic as part of Shovill assembly
                      (default: ${params.shovill_trim})
    --megahit_preset  Megahit paramaters preset. Possible values:
                      meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'
                      meta-large: '--k-min 27 --k-max 127 --k-step 10'
                      (large & complex metagenomes, like soil)

  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}

  ${c_bul}Other Options:${c_reset}
    ${c_green}--outdir${c_reset}          The output directory where the results will be saved
                      (default: ${c_green}${params.outdir}${c_reset})
    -w/--work-dir     The temporary directory where intermediate data will be 
                      saved (default: ${workflow.workDir})
    -profile          Configuration profile to use. [standard, singularity, 
                      conda, slurm] (default '${workflow.profile}')
    --tracedir        Pipeline run info output directory (default: 
                      ${params.tracedir})
  """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

//=============================================================================
// HELPER FUNCTIONS
//=============================================================================

def check_centrifuge_db(centrifuge_db) {
  file_centrifuge_db = file(centrifuge_db)
  prefix = file_centrifuge_db.getName()
  centrifuge_dir = file_centrifuge_db.getParent()
  if ( !centrifuge_dir.isDirectory() || !centrifuge_dir.exists() ) {
    exit 1, "Centrifuge DB does not exist at '$centrifuge_dir'! Please specify a valid Centrifuge DB."
  }
  any_valid = false
  centrifuge_dir.eachFile { f ->
    println "CENTRIFUGE FILE: $f"
    if ( f.isFile() ) {
      if ( f.getName().startsWith(prefix) && f.getExtension() == 'cf') {
        any_valid = true
      }
    }
  }
  if ( !any_valid ) {
    exit 1, "No valid Centrifuge DB files with prefix '$prefix' in '$centrifuge_dir' and extension 'cf'! Please specify a valid Centrifuge classification DB directory and prefix."
  }
}

def check_kraken2_db(kraken2_db) {
  kraken2_db_dir = file(kraken2_db)
  if ( !kraken2_db_dir.isDirectory() ) {
    exit 1, "The Kraken2 DB must be a directory! '$kraken2_db' is not a directory!"
  }
  if ( !kraken2_db_dir.exists() ) {
    exit 1, "The Kraken2 DB must be an existing directory! '$kraken2_db' does not exist!"
  }
}

def check_blastn_db(blastn_db) {
  if (blastn_db instanceof Boolean) {
    log.warn "--blastn_db not specified! Not running nucleotide BLAST on contigs!"
    return null
  }
  file_blastn_db = file(blastn_db)
  prefix = file_blastn_db.getName()
  db_dir = file_blastn_db.getParent()
  if ( !db_dir.isDirectory() || !db_dir.exists() ) {
    exit 1, "BLASTN DB directory does not exist at '${db_dir}'! Please specify a valid path to a BLAST DB."
  }
  any_valid = false
  nucl_db_extensions = ['nni', 'nhr', 'nin', 'nnd', 'nog', 'nsq'] as Set
  db_dir.eachFile { f -> 
    if ( f.isFile() ) {
      if ( f.getName().startsWith(prefix) && f.getExtension() in nucl_db_extensions) {
        any_valid = true
      }
    }
  }
  if ( !any_valid ) {
    exit 1, "Could not find valid nucleotide BLAST DB at '$blastn_db'. Please verify that this DB has files with extensions: $nucl_db_extensions"
  }
  return blastn_db
}


def check_megahit_preset(megahit_preset) {
  valid_presets = ['meta-sensitive', 'meta-large']
  if (megahit_preset instanceof Boolean) {
    return ''
  }
  if (megahit_preset instanceof String) {
    if (megahit_preset in valid_presets) {
      return megahit_preset
    }
  }
  log.warn "Invalid megahit_preset specified '${megahit_preset}'. Must be one of ${valid_presets}"
  return ''
}

// Check that all taxids are integers delimited by commas
def check_taxids(taxids) {
  if (taxids instanceof Boolean || taxids.toString().isEmpty()) {
    return null
  } else if (taxids.toString().isInteger()) {
    return taxids.toString()
  } else {
    taxids_list = taxids.toString()
      .split(',')
      .collect { it.strip() }
      .findAll { it != '' }
    if (!taxids_list.every { it.isInteger() }) {
      exit 1, "Not every element in `--taxids` is an integer!"
    }
    return taxids_list.join(',')
  }
}

// Get taxidlist for user-specified taxids
def taxidlist(taxids) {
  if (taxids instanceof String && file(taxids).exists()) {
    return file(taxids)
  }
  emptyFile = "EMPTY"
  if (taxids instanceof Boolean || taxids.toString().isEmpty()) {
    return emptyFile
  }
  l = taxids.toString()
                .split(',')
                .collect { it.strip() }
                .findAll { it != '' }
  if (l.size() == 0) {
    return emptyFile
  }
  if (l.size() == 1) {
    taxid = l[0].toInteger()
    if (taxid == 10239) {
      return file("$baseDir/data/Viruses-10239.taxidlist")
    } else if (taxid == 2) {
      return file("$baseDir/data/Bacteria-2.taxidlist")
    } else {
      return emptyFile
    }
  } else {
    return emptyFile
  }
}

//=============================================================================
// INPUT VALIDATION
//=============================================================================

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  exit 1, "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
}

taxids = check_taxids(params.taxids)
blastn_taxidlist = taxidlist(params.blastn_taxids)
megahit_preset = check_megahit_preset(params.megahit_preset)
if (params.centrifuge_db) check_centrifuge_db(params.centrifuge_db)
if (params.kraken2_db) check_kraken2_db(params.kraken2_db)
blastn_db = (params.blastn_db != null) ? check_blastn_db(params.blastn_db) : null



//=============================================================================
// WORKFLOW RUN PARAMETERS LOGGING
//=============================================================================

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']    = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Centrifuge DB'] = params.centrifuge_db ? params.centrifuge_db : "!NOT RUNNING CENTRIFUGE!"
summary['Kraken2 DB'] = params.kraken2_db ? params.kraken2_db : "!NOT RUNNING KRAKEN2!"
if(taxids && params.centrifuge_db && params.kraken2_db) {
  summary['Filter for reads in taxids'] = taxids
} else {
  summary['Filter reads by taxid?'] = 'NO'
}
summary['Output unclassified reads?'] = !(params.exclude_unclassified_reads as Boolean)
summary['Unicycler Mode']   = params.unicycler_mode
summary['Shovill Trim?']    = params.shovill_trim as Boolean
summary['Megahit Preset']   = megahit_preset
summary['BLASTN DB']        = blastn_db
summary['BLASTN taxidlist'] = blastn_taxidlist
if(params.blastn_qcov_hsp_perc != null) summary['BLASTN Q HSP %Cov'] = params.blastn_qcov_hsp_perc 
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = file(params.outdir)
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'slurm') summary['SLURM Queue'] = params.slurm_queue
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


if (params.container == null && (workflow.container ==~ /.*null/)) {
  exit 1, "Container null!"
}

//=============================================================================
// PROCESSES
//=============================================================================

// Remove any Coliphage phi-X174 reads using bbduk
process REMOVE_PHIX {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/phix_removal", pattern: "*.txt", mode: 'copy'
  publishDir "${params.outdir}/reads/phix_removed", pattern: "*.fastq.gz"
  
  input:
    file phix
    tuple val(sample_id), path(reads1), path(reads2)
  output:
    tuple val(sample_id), path(reads_out1), path(reads_out2), emit: 'reads'
    tuple val(sample_id), path(stats), emit: 'stats'

  script:
  reads_out1 = "${sample_id}_1.phix_removed.fastq.gz"
  reads_out2 = "${sample_id}_2.phix_removed.fastq.gz"
  stats = "${sample_id}-remove_phix-stats.txt"
  """
  bbduk.sh \\
    -Xmx${task.memory.toMega()}m \\
    in1=$reads1 in2=$reads2 \\
    out1=$reads_out1 out2=$reads_out2 \\
    ref=$phix k=31 hdist=1 \\
    stats=$stats
  """
}

// fastp adapter trimming and quality filtering
// By default, only up to 40% of bases can be below the min Phred score base
// quality threshold of Q15.
process FASTP {
  tag "$sample_id"
  publishDir "${params.outdir}/fastp/html", pattern: "*.html", mode: 'copy'
  publishDir "${params.outdir}/fastp/json", pattern: "*.json", mode: 'copy'
  publishDir "${params.outdir}/reads/fastp", pattern: "*.fastp.fastq.gz"

  input:
    tuple val(sample_id), path(r1), path(r2)
  output:
    tuple val(sample_id), path(reads_out1), path(reads_out2), emit: 'reads'
    tuple val(sample_id), path(html_report), path(json_report), emit: 'report'

  script:
  reads_out1 = "${sample_id}_1.fastp.fastq.gz"
  reads_out2 = "${sample_id}_2.fastp.fastq.gz"
  json_report = "fastp-report-${sample_id}.json"
  html_report = "fastp-report-${sample_id}.html"
  """
  fastp -i $r1 -I $r2 \\
    -o $reads_out1 -O $reads_out2 \\
    -p -c -R "$sample_id fastp report" \\
    -w ${task.cpus} \\
    -q ${params.fastp_min_base_quality} \\
    -u ${params.fastp_max_percent_low_qual_base} \\
    -j $json_report -h $html_report
  """
}

// FastQC
process FASTQC {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/fastqc", mode: 'copy',
      saveAs: { filename -> 
        filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
      }

  input:
    tuple val(sample_id), path(reads1), path(reads2)
  output:
    file "*_fastqc.{zip,html}"

  script:
  """
  fastqc -q $reads1 $reads2
  """
}

process KRAKEN2 {
  tag "$sample_id"
  publishDir "${params.outdir}/kraken2/results", pattern: "*-kraken2_results.tsv", mode: 'copy'
  publishDir "${params.outdir}/kraken2/reports", pattern: "*-kraken2_report.tsv", mode: 'copy'

  input:
    path(kraken2_db_dir)
    tuple val(sample_id),
          path(reads1),
          path(reads2)
  output:
    tuple val(sample_id),
          path(reads1),
          path(reads2),
          path(results),
          path(report)

  script:
  results = "${sample_id}-kraken2_results.tsv"
  report = "${sample_id}-kraken2_report.tsv"
  """
  kraken2 --memory-mapping --threads ${task.cpus} \\
    --db ./${kraken2_db_dir}/ \\
    --output ${results} \\
    --report ${report} \\
    $reads1 $reads2
  """
}

process CENTRIFUGE {
  tag "$sample_id"
  publishDir "${params.outdir}/centrifuge/results", pattern: "*-centrifuge_results.tsv", mode: 'copy'
  publishDir "${params.outdir}/centrifuge/reports", pattern: "*-centrifuge_kreport.tsv", mode: 'copy'
  // Memory usage as 125% of sum size of all cf index files for index multiplied
  memory {
    file_sizes = file(centrifuge_db_dir).listFiles()
      .findAll { it.isFile() && file(it).getExtension() == 'cf' }
      .collect { it.size() }
      .inject(0, { r, i -> r + i })
      .toLong()
    (file_sizes * 1.25).toLong()
  }

  input:
    tuple val(db_name), 
          path(centrifuge_db_dir)
    tuple val(sample_id),
          path(reads1),
          path(reads2)
  output:
    tuple val(sample_id),
          path(reads1),
          path(reads2),
          path(results),
          path(kreport)

  script:
  results = "${sample_id}-centrifuge_results.tsv"
  kreport = "${sample_id}-centrifuge_kreport.tsv"
  """
  centrifuge -x ${centrifuge_db_dir}/${db_name} \\
    -1 $reads1 -2 $reads2 \\
    -S $results -p ${task.cpus} 
  centrifuge-kreport -x ${centrifuge_db_dir}/${db_name} $results > $kreport
  """
}

process FILTER_READS_BY_CLASSIFICATIONS {
  tag "$sample_id"
  publishDir "${params.outdir}/filtered_reads/", pattern: "*.filtered.fastq.gz", mode: 'copy'

  input:
    tuple val(sample_id),
          path(reads1),
          path(reads2),
          path(kraken2_results),
          path(kraken2_report),
          path(centrifuge_results),
          path(centrifuge_report)
  output:
    tuple val(sample_id),
          path(filtered_reads1),
          path(filtered_reads2) optional true

  script:
  filtered_reads1 = "${sample_id}_1.filtered.fastq.gz"
  filtered_reads2 = "${sample_id}_2.filtered.fastq.gz"
  exclude_unclassified_reads = (params.exclude_unclassified_reads) ? '--exclude-unclassified' : ''
  """
  filter_classified_reads -i $reads1 -I $reads2 \\
    -o $filtered_reads1 -O $filtered_reads2 \\
    -c $centrifuge_results -C $centrifuge_report \\
    -k $kraken2_results -K $kraken2_report \\
    --taxids ${params.taxids} $exclude_unclassified_reads
  """
}

process UNICYCLER_ASSEMBLY {
  tag "$sample_id"
  publishDir "${params.outdir}/assemblies/unicycler/$sample_id", mode: 'copy'

  input:
    tuple(val(sample_id), path(reads1), path(reads2))
  output:
    tuple val(sample_id), val('unicycler'), path(output_contigs, optional: true), emit: 'contigs'
    tuple val(sample_id), path(output_unicycler_log, optional: true), emit: 'log'
    tuple val(sample_id), path(output_gfa, optional: true), emit: 'gfa'

  script:
  output_contigs = "${sample_id}-assembly.fasta"
  output_gfa = "${sample_id}-assembly.gfa"
  output_unicycler_log = "${sample_id}-unicycler.log"
  """
  unicycler -t ${task.cpus} --mode ${params.unicycler_mode} -o $sample_id -1 $reads1 -2 $reads2
  ln -s ${sample_id}/assembly.fasta $output_contigs
  ln -s ${sample_id}/assembly.gfa $output_gfa
  ln -s ${sample_id}/unicycler.log $output_unicycler_log
  """
}

process SHOVILL_ASSEMBLY {
  tag "$sample_id"
  publishDir "${params.outdir}/assemblies/shovill/$sample_id", mode: 'copy'

  input:
    tuple val(sample_id), path(reads1), path(reads2)
  output:
    tuple val(sample_id), val('shovill'), path(output_contigs, optional: true), emit: 'contigs'
    tuple val(sample_id), path(output_shovill_log, optional: true), emit: 'log'
    tuple val(sample_id), path(output_gfa, optional: true), emit: 'gfa'

  script:
  output_contigs = "${sample_id}-contigs.fasta"
  output_gfa = "${sample_id}-contigs.gfa"
  output_shovill_log = "${sample_id}-shovill.log"
  shovill_trim = (params.shovill_trim) ? "--trim" : ""
  """
  shovill --cpus ${task.cpus} --ram ${task.memory.toGiga()} \\
    --R1 $reads1 --R2 $reads2 \\
    --mincov 0.1 --depth 0 $shovill_trim \\
    --outdir $sample_id 
  ln -s ${sample_id}/contigs.fa $output_contigs
  ln -s ${sample_id}/contigs.gfa $output_gfa
  ln -s ${sample_id}/shovill.log $output_shovill_log
  """
}

process MEGAHIT_ASSEMBLY {
  publishDir "${params.outdir}/assemblies/megahit/$sample_id", mode: 'copy'

  input:
    tuple val(sample_id), path(reads1), path(reads2)
  output:
    tuple val(sample_id), val('megahit'), path(output_contigs, optional: true), emit: 'contigs'
    tuple val(sample_id), path(output_log, optional: true), emit: 'log'

  script:
  output_contigs = "${sample_id}-contigs.fasta"
  output_log = "${sample_id}-megahit.log"
  megahit_preset = (params.megahit_preset.toString() == '') ? '' : "--presets ${params.megahit_preset}"
  """
  megahit \\
    -t ${task.cpus} \\
    -m ${task.memory.toBytes()} \\
    $megahit_preset \\
    -1 $reads1 \\
    -2 $reads2 \\
    -o out \\
    --out-prefix out
  ln -s out/out.contigs.fa $output_contigs
  ln -s out/out.log $output_log
  """
}


process BLASTN {
  publishDir "${params.outdir}/blastn/$assembler", pattern: "*.tsv", mode: 'copy'
  tag "$sample_id|$dbname|$assembler|$txids"

  input:
    tuple val(dbname), path(dbdir)
    path(txids)
    tuple val(sample_id), val(assembler), path(contigs)
  output:
    tuple val(sample_id), val(assembler), path(contigs), path(blast_out)

  script:
  taxidlist_opt = (taxids == 'EMPTY') ? '' : "-taxidlist $txids"
  blastn_qcov_hsp_perc = (params.blastn_qcov_hsp_perc == null) ? '' : '-qcov_hsp_perc ${params.blastn_qcov_hsp_perc}'
  blast_out = "blastn-${sample_id}-VS-${dbname}.tsv"
  blast_tab_columns = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid ssciname"
  """
  blastn $taxidlist_opt $blastn_qcov_hsp_perc \\
    -num_threads ${task.cpus} \\
    -db $dbdir/$dbname \\
    -query $contigs \\
    -outfmt "6 $blast_tab_columns" \\
    -evalue 1e-6 \\
    -out $blast_out
  """
}


process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path('*') 
    path(config) 

    output:
    path('multiqc_report.html')

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc .
    """
}

//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================
workflow {
  // Create channel for paired end reads
  ch_reads = Channel.fromFilePairs(
      params.reads,
      flat: true,
      checkIfExists: true)
    .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}*.fastq.gz\"' (quotes around reads path required if using `*` and other characters expanded by the shell!)"}
    .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1], it[2] ] }

  // Value channel with file containing phix sequence from workflow data folder
  ch_phix = Channel.value(file("$baseDir/data/phix.fa"))
  REMOVE_PHIX(ch_phix, ch_reads)
  ch_qced_reads = REMOVE_PHIX.out.reads
  if (!params.skip_fastp) {
    FASTP(REMOVE_PHIX.out.reads)
    ch_qced_reads = FASTP.out.reads
  }
  fastqc_reports = FASTQC(ch_qced_reads)

  if (params.kraken2_db) {
    KRAKEN2(Channel.value(file(params.kraken2_db)), ch_qced_reads)
  }
  if (params.centrifuge_db) {
    cdb = file(params.centrifuge_db)
    ch_centrifuge_db = Channel.value([cdb.getName(), cdb.getParent()])
    CENTRIFUGE(ch_centrifuge_db, ch_qced_reads)
  }
  if (params.kraken2_db && params.centrifuge_db && taxids) {
    ch_kraken2_and_centrifuge_results = KRAKEN2.out
      .join(CENTRIFUGE.out, remainder: true)
      .map { sample_id, r1, r2, kraken2_results, kraken2_report, _r1, _r2, centrifuge_results, centrifuge_kreport -> 
        [sample_id, r1, r2, kraken2_results, kraken2_report, centrifuge_results, centrifuge_kreport]
      }
    ch_kraken2_and_centrifuge_results \
      | FILTER_READS_BY_CLASSIFICATIONS \
      | (UNICYCLER_ASSEMBLY & SHOVILL_ASSEMBLY & MEGAHIT_ASSEMBLY)
  } else {
    ch_qced_reads | (UNICYCLER_ASSEMBLY & SHOVILL_ASSEMBLY & MEGAHIT_ASSEMBLY)
  }

  // Run BLASTN if valid BLAST DB specified
  if (blastn_db != null) {
    // BLASTN input channels
    file_blastn_db = file(blastn_db)
    ch_blastdb = Channel.value( [file_blastn_db.getName(), file_blastn_db.getParent()] )
    ch_taxidlist = Channel.value( blastn_taxidlist )
    ch_contigs = UNICYCLER_ASSEMBLY.out.contigs
      .mix(SHOVILL_ASSEMBLY.out.contigs, MEGAHIT_ASSEMBLY.out.contigs)

    BLASTN(ch_blastdb, ch_taxidlist, ch_contigs)
  }
}

//=============================================================================
// WORKFLOW EVENT HANDLERS
//=============================================================================
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
