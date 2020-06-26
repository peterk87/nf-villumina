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
