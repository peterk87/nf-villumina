process BLASTN {
  publishDir "${params.outdir}/blastn/$assembler", pattern: "*.tsv", mode: 'copy'
  tag "$sample_id|$dbname|$assembler|$txids"

  input:
    tuple dbname, path(dbdir)
    path(txids)
    tuple val(sample_id), val(assembler), path(contigs)
  output:
    tuple val(sample_id), val(assembler), path(contigs), path(blast_out)

  script:
  taxidlist_opt = (taxids == 'EMPTY') ? '' : "-taxidlist $txids"
  blast_out = "blastn-${sample_id}-VS-${dbname}.tsv"
  blast_tab_columns = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid ssciname"
  """
  blastn $taxidlist_opt \\
    -num_threads ${task.cpus} \\
    -db $dbdir/$dbname \\
    -query $contigs \\
    -outfmt "6 $blast_tab_columns" \\
    -evalue 1e-6 \\
    -out $blast_out
  """
}
