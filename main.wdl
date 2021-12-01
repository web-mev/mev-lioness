workflow mevLioness {
    
    # A fixed motif file
    File motif_file
    # A fixed PPI file
    File ppi_file
    # A user uploaded exprs count matrix
    File exprs_file

    # Number of scatters to split the samples over
    Int num_ranges = 10

    call runPandas {
        input:
            motif_file = motif_file,
            ppi_file = ppi_file,
            exprs_file = exprs_file,
            num_ranges = num_ranges
    }

    scatter (line_num in range(num_ranges)) {
        call runLioness {
            input:
                panda_pickle = runPandas.panda_pickle,
                slice_tsv = runPandas.scatter_ranges,
                line_num = line_num
        }
    }

}

task runPandas {
    File motif_file
    File ppi_file
    File exprs_file

    Int disk_size = 40

    command {
        python3 /opt/software/panda.py \
            --motif ${motif_file} \
            --ppi ${ppi_file} \
            ${exprs_file}
    }

    output {
        File panda_pickle = "panda_obj.pkl"
        File panda_matrix = "panda_output.mtx"
        File scatter_ranges = "sample_scatter_ranges.txt"
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
            ${panda_pickle}
    }

    output {
        File lioness_matrix = "lioness_matrix.npy"
    }

    runtime {
        docker: "hsphqbrc/mev-netzoopy"
        cpu: 8
        memory: "128 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}