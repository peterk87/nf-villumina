
process CENTRIFUGE_INSPECT_INDEX {
  tag "$centrifuge_db_dir / $db_name"

  input:
  tuple db_name, path(centrifuge_db_dir)
  output:
  tuple path(name_table), path(taxonomy_tree)

  script:
  name_table = "centrifuge-${db_name}-name_table"
  taxonomy_tree = "centrifuge-${db_name}-taxonomy_tree"
  """
  centrifuge-inspect --name-table ${centrifuge_db_dir}/${db_name} > $name_table
  centrifuge-inspect --taxonomy-tree ${centrifuge_db_dir}/${db_name} > $taxonomy_tree
  """
}

process CENTRIFUGE {
  tag "$sample_id"
  publishDir "${params.outdir}/centrifuge/results", pattern: "*-centrifuge_results.tsv", mode: 'copy'
  // Memory usage as 125% of sum size of all cf index files for index multiplied
  memory {
    file_sizes = file(centrifuge_db_dir).listFiles()
      .findAll { it.isFile() && file(it).getExtension() == 'cf' }
      .collect { it.size() }
      .inject(0, { r, i -> r + i })
      .toLong()
    (file_sizes * 1.25).toLong()
  }

  input:
  tuple db_name, 
        path(centrifuge_db_dir)
  tuple sample_id,
        path(reads1),
        path(reads2)
  output:
  tuple sample_id,
        path(reads1),
        path(reads2),
        path(results)

  script:
  results = "${sample_id}-centrifuge_results.tsv"
  kreport = "${sample_id}-centrifuge_kreport.tsv"
  """
  centrifuge -x ${centrifuge_db_dir}/${db_name} \\
    -1 $reads1 -2 $reads2 \\
    -S $results -p ${task.cpus}
  """
}

process CENTRIFUGE_KREPORT {
  tag "$sample_id"
  publishDir "${params.outdir}/centrifuge/reports", pattern: "*-centrifuge_kreport.tsv", mode: 'copy'

  input:
  tuple sample_id,
      path(reads1),
      path(reads2),
      path(results)
  tuple path(name_table), path(taxonomy_tree)
  
  output:
  tuple sample_id,
      path(reads1),
      path(reads2),
      path(results),
      path(kreport)

  script:
  """
  centrifuge-kreport-mod.pl --names $results > $kreport
  """
}
