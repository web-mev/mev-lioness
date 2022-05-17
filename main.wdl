workflow mevLioness {
    
    # A fixed output name for scatter ranges
    String scatter_ranges_fname = "sample_scatter_ranges.txt"

    # Maximum number of samples allowed in a slice
    Int max_num_in_slice

    # A fixed motif file
    File motif_file
    # A fixed PPI file
    File ppi_file
    # A user uploaded exprs count matrix
    File exprs_file

    # Output PANDA matrix name
    String panda_matrix_filename


    call determineScatters {
        input: 
            exprs_file = exprs_file,
            max_num_in_slice = max_num_in_slice
    }

    call runPandas {
        input:
            scatter_ranges_fname = scatter_ranges_fname,
            num_ranges = determineScatters.num_ranges,
            motif_file = motif_file,
            ppi_file = ppi_file,
            panda_matrix_filename = panda_matrix_filename,
            exprs_file = exprs_file
    }

    scatter (line_num in range(determineScatters.num_ranges)) {
        call runLioness {
            input:
                panda_pickle = runPandas.panda_pickle,
                slice_tsv = runPandas.scatter_ranges,
                line_num = line_num,
                exprs_file = exprs_file
        }
    }

    call mergeLioness {
        input:
            unrolled_output_filename = unrolled_output_filename,
            gene_ts_output_filename = gene_ts_output_filename,
            tf_ts_output_filename = tf_ts_output_filename,
            lioness_scatter_tsv = runLioness.lioness_scatter_tsv
    }

    output {
        File panda_tsv = runPandas.panda_matrix
        File lioness_gene_ts_tsv = mergeLioness.lioness_gene_target_scores
        File lioness_tf_ts_tsv = merge_lioness.lioness_tf_target_scores
        File lioness_unrolled_tsv = merge_lioness.lioness_unrolled
    }
}

task determineScatters {
    File exprs_file
    Int max_num_in_slice
    
    command {
        python3 /opt/software/determine_scatter.py \
            --max ${max_num_in_slice} \
            ${exprs_file}
    }

    output {
        Int num_ranges = stdout()
    }
}

task runPandas {
    String scatter_ranges_fname
    Int num_ranges
    File motif_file
    File ppi_file
    String panda_matrix_filename
    File exprs_file

    Int disk_size = 40

    command {
        python3 /opt/software/panda.py \
            --scatter ${scatter_ranges_fname} \
            --ranges ${num_ranges} \
            --motif ${motif_file} \
            --ppi ${ppi_file} \
            --out ${panda_matrix_filename} \
            ${exprs_file}
    }

    output {
        File panda_pickle = "panda_obj.pkl"
        File panda_matrix = "${panda_matrix_filename}"
        File scatter_ranges = "${scatter_ranges_fname}"
    }

    runtime {
        docker: "hsphqbrc/mev-netzoopy"
        cpu: 8
        memory: "128 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task runLioness {
    File panda_pickle
    File slice_tsv
    Int line_num
    
    Int disk_size = 40

    command {
        python3 /opt/software/lioness.py \
            --slices ${slice_tsv} \
            --line ${line_num} \
            --exprs ${exprs_file} \
            ${panda_pickle}
    }

    output {
        File lioness_scatter_tsv = "lioness_scatter_output.tsv"
    }

    runtime {
        docker: "hsphqbrc/mev-netzoopy"
        cpu: 8
        memory: "128 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task mergeLioness {
    String unrolled_output_filename
    String gene_ts_output_filename
    String tf_ts_output_filename
    Array[File] lioness_scatter_tsv

    Int disk_size = 40

    command {
        python3 /opt/software/merge_lioness.py \
            --full ${unrolled_output_filename} \
            --gene ${gene_ts_output_filename} \
            --tf ${tf_ts_output_filename} \
            --lioness ${sep=" " lioness_scatter_tsv};
    }

    output {
        File lioness_gene_target_scores = "${gene_ts_output_filename}"
        File lioness_tf_target_scores = "${tf_ts_output_filename}"
        File lioness_unrolled = "${unrolled_output_filename}"
    }

    runtime {
        docker: "hsphqbrc/mev-netzoopy"
        cpu: 4
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
