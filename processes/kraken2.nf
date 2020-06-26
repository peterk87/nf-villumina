process KRAKEN2 {
  tag "$sample_id"
  publishDir "${params.outdir}/kraken2/results", pattern: "*-kraken2_results.tsv", mode: 'copy'
  publishDir "${params.outdir}/kraken2/reports", pattern: "*-kraken2_report.tsv", mode: 'copy'

  input:
    path(kraken2_db_dir)
    tuple sample_id,
          path(reads1),
          path(reads2)
  output:
    tuple sample_id,
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
