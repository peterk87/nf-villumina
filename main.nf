#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ==================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ==================================================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run --reads ${params.reads} --outdir ${params.outdir}

    Options:
      --reads                       Input reads directory and pattern (default: "${params.reads}")
      --centrifuge_db               Path to Centrifuge DB and prefix (default: ${params.centrifuge_db})
      --kraken2_db                  Path to Kraken2 DB (default: ${params.kraken2_db})
      --unicycler_mode              Unicycler assembly mode (default: ${params.unicycler_mode})
    Other options:
      --outdir                      The output directory where the results will be saved (default: ${params.outdir})
      -w/--work-dir                 The temporary directory where intermediate data will be saved (default: ${workflow.workDir})
      -profile                      Configuration profile to use. [standard, other_profiles] (default '${workflow.profile}')
      --tracedir                    Pipeline run info output directory (default: ${params.tracedir})
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

def check_centrifuge_db(centrifuge_db) {
  file_centrifuge_db = file(centrifuge_db)
  prefix = file_centrifuge_db.getName()
  centrifuge_dir = file_centrifuge_db.getParent()
  if ( !centrifuge_dir.isDirectory() || !centrifuge_dir.exists() ) {
    exit 1, "Centrifuge DB does not exist at '$centrifuge_dir'! Please specify a valid Centrifuge DB."
  }
  any_valid = false
  centrifuge_dir.eachFile { f ->
    if ( f.isFile() ) {
      if ( f.getName() =~ /^$prefix/ && f.getExtension() == 'cf') {
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

if (workflow.containerEngine == 'singularity') {

}

check_centrifuge_db(params.centrifuge_db)
check_kraken2_db(params.kraken2_db)

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
summary['Pipeline Name']  = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Centrifuge DB'] = params.centrifuge_db
summary['Kraken2 DB']   = params.kraken2_db
summary['Unicycler Mode'] = params.unicycler_mode
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

outdir = params.outdir

Channel.value( 
  [ 
    file(params.centrifuge_db).getName(), 
    file(params.centrifuge_db).getParent() 
  ] )
  .set { ch_centrifuge_db }

Channel.value( file(params.kraken2_db) )
  .set { ch_kraken2_db }

Channel
  .fromFilePairs( params.reads, flat: true)
  .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '~/my-reads/*R{1,2}*.fastq.gz'"}
  .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1], it[2] ] }
  .dump(tag: 'ch_reads')
  .set { ch_reads }


process download_phix {
  output:
    file('phix.fa') into phix_file

  '''
  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text" > phix.fa
  '''
}

process remove_phix {
  tag "$sample_id"
  
  input:
    file phix from phix_file
    set sample_id, file(read1), file(read2) from ch_reads
  output:
    set sample_id, file("R1.fq"), file("R2.fq") into ch_reads_minus_phix
    set sample_id, file("${sample_id}-remove_phix-stats.txt") into ch_remove_phix_stats

  script:
  """
  bbduk.sh in1=${read1} in2=${read2} out1=R1.fq out2=R2.fq ref=$phix k=31 hdist=1 stats=${sample_id}-remove_phix-stats.txt
  """
}

process fastp {
  tag "$sample_id"
  publishDir "$outdir/fastp/html", pattern: "*.html", mode: 'copy'
  publishDir "$outdir/fastp/json", pattern: "*.json", mode: 'copy'

  input:
    set sample_id, file(r1), file(r2) from ch_reads_minus_phix
  output:
    set sample_id, file(reads_out1), file(reads_out2) into ch_reads_for_fastqc, ch_reads_for_centrifuge, ch_reads_for_kraken2
    set sample_id, file(html_report), file(json_report)

  script:
  reads_out1 = "${sample_id}_1.fastq"
  reads_out2 = "${sample_id}_2.fastq"
  json_report = "fastp-report-${sample_id}.json"
  html_report = "fastp-report-${sample_id}.html"
  """
  fastp -i $r1 -I $r2 -o $reads_out1 -O $reads_out2 -p -c -R "$sample_id fastp report" -w ${task.cpus} -q 15 -j $json_report -h $html_report
  """
}

/*
process index_fastqs {
  tag "$sample_id"

  input:
    set sample_id, file(r1), file(r2) from ch_reads_for_faidx
  output:
    set sample_id, file(r1), file(r2), file('*.fai') into ch_faidx

  script:
  """
  samtools faidx $r1
  samtools faidx $r2
  """
}
*/

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(sample_id), file(reads1), file(reads2) from ch_reads_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q $reads1 $reads2
    """
}

// TODO: configure MultiQC
/*
 * STEP 2 - MultiQC
 */
/*
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
      //file multiqc_config from ch_multiqc_config
      // TODO nf-core: Add in log files from your new processes for MultiQC to find!
      file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
      // file ('software_versions/*') from software_versions_yaml
      //file workflow_summary from create_workflow_summary(summary)

    output:
      file "*multiqc_report.html" into multiqc_report
      file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename .
    """
}
*/

process kraken2_classification {
  tag "$sample_id"
  publishDir "$outdir/kraken2/results", pattern: "*-kraken2_results.tsv", mode: 'copy'
  publishDir "$outdir/kraken2/reports", pattern: "*-kraken2_report.tsv", mode: 'copy'

  input:
    file(kraken2_db_dir) from ch_kraken2_db
    set val(sample_id), file(reads1), file(reads2) from ch_reads_for_kraken2
  output:
    set val(sample_id), file(reads1), file(reads2), file(results), file(report) into ch_kraken2

  script:
  results = "${sample_id}-kraken2_results.tsv"
  report = "${sample_id}-kraken2_report.tsv"
  """
  kraken2 --memory-mapping --db ${kraken2_db_dir} --threads ${task.cpus} --output ${results} --report ${report} $reads1 $reads2
  """
}

process centrifuge {
  tag "$sample_id"
  publishDir "$outdir/centrifuge/$sample_id", pattern: "*.tsv", mode: 'copy'

  input:
    set db_name, file(centrifuge_db_dir) from ch_centrifuge_db
    set val(sample_id), file(reads1), file(reads2) from ch_reads_for_centrifuge
  output:
    set val(sample_id), file(reads1), file(reads2), file(results), file(report) into ch_centrifuge_results

  script:
  results = "${sample_id}-centrifuge_results.tsv"
  report = "${sample_id}-report.tsv"
  """
  centrifuge -x ${centrifuge_db_dir}/${db_name} -1 $reads1 -2 $reads2 -S $results --report-file $report --mm
  """
}

process centrifuge_kraken_report {
  tag "$sample_id"
  publishDir "$outdir/centrifuge", pattern: "*-kreport.tsv", mode: 'copy'

  input:
    set db_name, file(centrifuge_db_dir) from ch_centrifuge_db
    set val(sample_id), file(reads1), file(reads2), file(results), file(report) from ch_centrifuge_results
  output:
    set val(sample_id), file(reads1), file(reads2), file(results), file(kreport) into ch_centrifuge_kraken_report_results

  script:
  kreport = "${sample_id}-kreport.tsv"
  """
  centrifuge-kreport -x ${centrifuge_db_dir}/${db_name} $results > $kreport
  """
}

ch_kraken2
  .join(ch_centrifuge_kraken_report_results, remainder: true)
  .map { sample_id, reads1, reads2, kraken2_results, kraken2_report, r1, r2, centrifuge_results, centrifuge_kreport -> 
    [sample_id, reads1, reads2, kraken2_results, kraken2_report, centrifuge_results, centrifuge_kreport]
  }
  .dump(tag: "ch_kraken2_and_centrifuge")
  .set { ch_kraken2_and_centrifuge }

process filter_reads_by_classifications {
  tag "$sample_id"
  publishDir "$outdir/filtered_reads/", pattern: "*.viral_unclassified.fastq", mode: 'copy'

  input:
    set sample_id, file(reads1), file(reads2), file(kraken2_results), file(kraken2_report), file(centrifuge_results), file(centrifuge_report) from ch_kraken2_and_centrifuge
  output:
    set sample_id, file(filtered_reads1), file(filtered_reads2) optional true into ch_reads_for_unicycler, ch_reads_for_shovill

  script:
  filtered_reads1 = "${sample_id}_1.viral_unclassified.fastq"
  filtered_reads2 = "${sample_id}_2.viral_unclassified.fastq"
  """
  filter_classified_reads -i $reads1 -I $reads2 -o $filtered_reads1 -O $filtered_reads2 -c $centrifuge_results -C $centrifuge_report -k $kraken2_results -K $kraken2_report
  """
}

process unicycler_assembly {
  tag "$sample_id"
  publishDir "$outdir/assemblies/unicycler/$sample_id", mode: 'copy'

  input:
    set val(sample_id), file(reads1), file(reads2) from ch_reads_for_unicycler
  output:
    set val(sample_id), val('unicycler'), file(output_contigs) optional true into ch_unicycler_assembly
    set sample_id, file(output_unicycler_log) optional true into ch_unicycler_assembly_log
    set sample_id, file(output_gfa) optional true into ch_unicycler_assembly_gfa

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

process shovill_assembly {
  tag "$sample_id"
  publishDir "$outdir/assemblies/shovill/$sample_id", mode: 'copy'

  input:
    set val(sample_id), file(reads1), file(reads2) from ch_reads_for_shovill
  output:
    set val(sample_id), val('shovill'), file(output_contigs) optional true into ch_shovill_assembly
    set sample_id, file(output_shovill_log) optional true into ch_shovill_assembly_log
    set sample_id, file(output_gfa) optional true into ch_shovill_assembly_gfa

  script:
  output_contigs = "${sample_id}-contigs.fasta"
  output_gfa = "${sample_id}-contigs.gfa"
  output_shovill_log = "${sample_id}-shovill.log"
  """
  shovill --R1 $reads1 --R2 $reads2 --cpus ${task.cpus} --mincov 0.1 --depth 0 --outdir $sample_id --trim
  ln -s ${sample_id}/contigs.fa $output_contigs
  ln -s ${sample_id}/contigs.gfa $output_gfa
  ln -s ${sample_id}/shovill.log $output_shovill_log
  """
}

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
