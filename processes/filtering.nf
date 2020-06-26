process FILTER_READS_BY_CLASSIFICATIONS {
  tag "$sample_id"
  publishDir "${params.outdir}/filtered_reads/", pattern: "*.filtered.fastq.gz", mode: 'copy'

  input:
    tuple sample_id,
          path(reads1),
          path(reads2),
          path(kraken2_results),
          path(kraken2_report),
          path(centrifuge_results),
          path(centrifuge_report)
  output:
    tuple sample_id,
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
