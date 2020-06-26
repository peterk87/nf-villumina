// Remove any Coliphage phi-X174 reads using bbduk
process REMOVE_PHIX {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/phix_removal", pattern: "*.txt", mode: 'copy'
  publishDir "${params.outdir}/reads/phix_removed", pattern: "*.fastq.gz"
  
  input:
    file phix
    tuple sample_id, path(reads1), path(reads2)
  output:
    tuple sample_id, path(reads_out1), path(reads_out2), emit: 'reads'
    tuple sample_id, path(stats), emit: 'stats'

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
