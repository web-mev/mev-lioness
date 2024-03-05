// a few hardcoded values:
params.max_num_in_slice = 10
params.motif_files_map = [
    symbol: "s3://webmev-public/tissues_motif.symbol.tsv",
    ensembl: "s3://webmev-public/tissues_motif.ensg.tsv"
]
params.ppi_file = "s3://webmev-public/tissues_ppi.tsv"

process run_lioness {

    tag "Run LIONESS on shard"
    container "ghcr.io/web-mev/mev-lioness"
    cpus 8
    memory '128 GB'

    input:
        path exprs_file
        path panda_pickle
        tuple val(start), val(end)

    output:
        path "${output_fname}"

    script:
        tmp_dir = "/tmp_data"
        output_fname = "lioness_scatter_output.${start}_${end}.tsv"
        """
        mkdir ${tmp_dir}

        python3 /usr/local/bin/lioness.py \
            --start ${start} \
            --end ${end} \
            --exprs ${exprs_file} \
            --output ${output_fname} \
            --save_dir ${tmp_dir} \
            ${panda_pickle}
        """
}

process determine_scatters {
    tag "LIONESS scatter"
    container "ghcr.io/web-mev/mev-lioness"
    cpus 2
    memory '8 GB'

    input:
        path exprs_file
        val max_num_in_slice

    output:
        stdout

    script:
        """
        python3 /usr/local/bin/determine_scatter.py \
            --max ${max_num_in_slice} \
            ${exprs_file}
        """
}

process run_panda {
    tag "Run panda"
    container "ghcr.io/web-mev/mev-lioness"
    cpus 8
    memory 120 GB

    input:
        path exprs_file
        val num_scatters

    output:
        path "${panda_pkl}"
        path "${scatter_ranges_fname}"

    script:
        ppi_file = params.ppi_file
        motif_file = params.motif_files_map["${params.identifier_choice}"]
        panda_pkl = "panda_obj.pkl"
        scatter_ranges_fname = "sample_scatter_ranges.csv"
        """
        python3 /usr/local/bin/panda.py \
            --scatter ${scatter_ranges_fname} \
            --num_scatters ${num_scatters} \
            --motif ${motif_file} \
            --ppi ${ppi_file} \
            ${exprs_file}
        """
}

process merge_lioness_shards {
    tag "Merge the LIONESS shards"
    publishDir "${params.output_dir}/mevLioness.lioness_gene_ts_tsv", mode:"copy", pattern:"${gene_ts_output_filename}"
    publishDir "${params.output_dir}/mevLioness.lioness_tf_ts_tsv", mode:"copy", pattern:"${tf_ts_output_filename}"
    container "ghcr.io/web-mev/mev-lioness"
    cpus 4
    memory 16 GB

    input:
        path 'shards'

    output:
        path "${gene_ts_output_filename}"
        path "${tf_ts_output_filename}"

    script:
        gene_ts_output_filename = "lioness_gene_target_scores.tsv"
        tf_ts_output_filename = "lioness_transcription_factor_target_scores.tsv"
        """
        python3 /usr/local/bin/merge_lioness.py \
            --gene ${gene_ts_output_filename} \
            --tf ${tf_ts_output_filename} \
            --lioness shards*
        """  
}

workflow {

    num_scatter_ch = determine_scatters(params.exprs_file, params.max_num_in_slice)
    (pkl_ch, scatter_range_ch) = run_panda(params.exprs_file, num_scatter_ch)
    ranges_ch = scatter_range_ch.splitCsv(header: true).map{
        row -> tuple(row.start, row.end)
    }
    // need to run 'collect' so that it will appropriately scatter the same
    // file across the lioness shards
    pkl_ch = pkl_ch.collect()

    // Run LIONESS shard-by-shard. The `collect` takes the channels and 
    // puts them together so that the merging takes all the files. Otherwise
    // you get a merge run on each individual file
    lioness_shards_ch = run_lioness(params.exprs_file, pkl_ch, ranges_ch).collect()

    merge_lioness_shards(lioness_shards_ch)
}
