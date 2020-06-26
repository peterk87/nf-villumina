#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include { 
  helpMessage; 
  check_taxids; 
  taxidlist; 
  check_megahit_preset; 
  check_centrifuge_db; 
  check_kraken2_db; 
  check_blastn_db 
  } from './lib/helpers'

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
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

include REMOVE_PHIX from './processes/remove_phix'
include FASTP from './processes/fastp'
include FASTQC from './processes/fastqc'
include KRAKEN2 from './processes/kraken2'
include { CENTRIFUGE_INSPECT_INDEX; CENTRIFUGE; CENTRIFUGE_KREPORT } from './processes/centrifuge'
include FILTER_READS_BY_CLASSIFICATIONS from './processes/filtering'
include { MEGAHIT_ASSEMBLY; SHOVILL_ASSEMBLY; UNICYCLER_ASSEMBLY } from './processes/assembly'
include BLASTN from './processes/blast'

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
  FASTP(REMOVE_PHIX.out.reads)
  fastqc_reports = FASTQC(FASTP.out.reads)

  if (params.kraken2_db) {
    KRAKEN2(Channel.value(file(params.kraken2_db)), FASTP.out.reads)
  }
  if (params.centrifuge_db) {
    cdb = file(params.centrifuge_db)
    ch_centrifuge_db = Channel.value([cdb.getName(), cdb.getParent()])
    CENTRIFUGE_INSPECT_INDEX(ch_centrifuge_db)
    CENTRIFUGE(ch_centrifuge_db, FASTP.out.reads)
    CENTRIFUGE_KREPORT(CENTRIFUGE.out, CENTRIFUGE_INSPECT_INDEX.out)
  }
  if (params.kraken2_db && params.centrifuge_db && taxids) {
    ch_kraken2_and_centrifuge_results = KRAKEN2.out
      .join(CENTRIFUGE_KREPORT.out, remainder: true)
      .map { sample_id, r1, r2, kraken2_results, kraken2_report, _r1, _r2, centrifuge_results, centrifuge_kreport -> 
        [sample_id, r1, r2, kraken2_results, kraken2_report, centrifuge_results, centrifuge_kreport]
      }
    ch_kraken2_and_centrifuge_results \
      | FILTER_READS_BY_CLASSIFICATIONS \
      | (UNICYCLER_ASSEMBLY & SHOVILL_ASSEMBLY & MEGAHIT_ASSEMBLY)
  } else {
    FASTP.out.reads | (UNICYCLER_ASSEMBLY & SHOVILL_ASSEMBLY & MEGAHIT_ASSEMBLY)
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
