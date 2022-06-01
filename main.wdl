workflow mevLioness {
    
    # A fixed output name for scatter ranges
    String scatter_ranges_fname = "sample_scatter_ranges.txt"

    # Maximum number of samples allowed in a slice
    Int max_num_in_slice = 10

    # A fixed motif file
    File? motif_file = "gs://mev-public-data/tissues_motif.tsv"

    # A fixed PPI file
    File? ppi_file = "gs://mev-public-data/tissues_ppi.tsv"
    
    # A name for the output file containing the target
    # scores for the genes
    String gene_ts_output_filename = "lioness_gene_target_scores.tsv"

    # A name for the output file containing the target
    # scores for the transcription factors
    String tf_ts_output_filename = "lioness_transcription_factor_target_scores.tsv"

    # A user uploaded exprs count matrix
    File exprs_file


    call determineScatters {
        input: 
            exprs_file = exprs_file,
            max_num_in_slice = max_num_in_slice
    }

    call runPanda {
        input:
            scatter_ranges_fname = scatter_ranges_fname,
            num_scatters = determineScatters.num_scatters,
            motif_file = motif_file,
            ppi_file = ppi_file,
            exprs_file = exprs_file
    }

    scatter (line_num in range(determineScatters.num_scatters)) {
        call runLioness {
            input:
                panda_pickle = runPanda.panda_pickle,
                slice_tsv = runPanda.scatter_ranges,
                line_num = line_num,
                exprs_file = exprs_file
        }
    }

    call mergeLioness {
        input:
            gene_ts_output_filename = gene_ts_output_filename,
            tf_ts_output_filename = tf_ts_output_filename,
            lioness_scatter_tsv = runLioness.lioness_scatter_tsv
    }

    output {
        File lioness_gene_ts_tsv = mergeLioness.lioness_gene_target_scores
        File lioness_tf_ts_tsv = mergeLioness.lioness_tf_target_scores
    }
}

task determineScatters {
    File exprs_file
    Int max_num_in_slice
    
    command {
        python3 /opt/software/determine_scatter.py \
            --max ${max_num_in_slice} \
            ${exprs_file} > nscatter.txt
    }

    output {
        Int num_scatters = read_int("nscatter.txt")
    }
}

task runPanda {
    String scatter_ranges_fname
    Int num_scatters
    File motif_file
    File ppi_file
    File exprs_file

    Int disk_size = 40

    command {
        python3 /opt/software/panda.py \
            --scatter ${scatter_ranges_fname} \
            --num_scatters ${num_scatters} \
            --motif ${motif_file} \
            --ppi ${ppi_file} \
            ${exprs_file}
    }

    output {
        File panda_pickle = "panda_obj.pkl"
        File scatter_ranges = "${scatter_ranges_fname}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-lioness"
        cpu: 8
        memory: "128 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task runLioness {
    File panda_pickle
    File exprs_file
    File slice_tsv
    Int line_num
    
    String output_fname = "lioness_scatter_output.tsv"
    String tmp_dir = "/tmp_data"
    Int disk_size = 40

    command {

        mkdir ${tmp_dir}

        python3 /opt/software/lioness.py \
            --slices ${slice_tsv} \
            --line ${line_num} \
            --exprs ${exprs_file} \
            --output ${output_fname} \
            --save_dir ${tmp_dir} \
            ${panda_pickle}
    }

    output {
        File lioness_scatter_tsv = "${output_fname}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-lioness"
        cpu: 8
        memory: "128 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task mergeLioness {
    String gene_ts_output_filename
    String tf_ts_output_filename
    Array[File] lioness_scatter_tsv

    Int disk_size = 40

    command {
        python3 /opt/software/merge_lioness.py \
            --gene ${gene_ts_output_filename} \
            --tf ${tf_ts_output_filename} \
            --lioness ${sep=" " lioness_scatter_tsv};
    }

    output {
        File lioness_gene_target_scores = "${gene_ts_output_filename}"
        File lioness_tf_target_scores = "${tf_ts_output_filename}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-lioness"
        cpu: 4
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
