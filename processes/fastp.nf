// fastp adapter trimming and quality filtering
// By default, only up to 40% of bases can be below the min Phred score base
// quality threshold of Q15.
process FASTP {
  tag "$sample_id"
  publishDir "${params.outdir}/fastp/html", pattern: "*.html", mode: 'copy'
  publishDir "${params.outdir}/fastp/json", pattern: "*.json", mode: 'copy'
  publishDir "${params.outdir}/reads/fastp", pattern: "*.fastp.fastq.gz"

  input:
    tuple sample_id, path(r1), path(r2)
  output:
    tuple sample_id, path(reads_out1), path(reads_out2), emit: 'reads'
    tuple sample_id, path(html_report), path(json_report), emit: 'report'

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
