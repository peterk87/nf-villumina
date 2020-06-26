process UNICYCLER_ASSEMBLY {
  tag "$sample_id"
  publishDir "${params.outdir}/assemblies/unicycler/$sample_id", mode: 'copy'

  input:
    tuple(val(sample_id), path(reads1), path(reads2))
  output:
    tuple sample_id, val('unicycler'), path(output_contigs, optional: true), emit: 'contigs'
    tuple sample_id, path(output_unicycler_log, optional: true), emit: 'log'
    tuple sample_id, path(output_gfa, optional: true), emit: 'gfa'

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
    tuple sample_id, val('shovill'), path(output_contigs, optional: true), emit: 'contigs'
    tuple sample_id, path(output_shovill_log, optional: true), emit: 'log'
    tuple sample_id, path(output_gfa, optional: true), emit: 'gfa'

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
    tuple sample_id, val('megahit'), path(output_contigs, optional: true), emit: 'contigs'
    tuple sample_id, path(output_log, optional: true), emit: 'log'

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
