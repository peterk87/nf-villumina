#!/usr/bin/env nextflow

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