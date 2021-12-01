workflow mevLioness {
    
    # A fixed motif file
    File motif_file
    # A fixed PPI file
    File ppi_file
    # A user uploaded exprs count matrix
    File exprs_file

    call runPandas {
        input:
            motif_file = motif_file,
            ppi_file = ppi_file,
            exprs_file = exprs_file
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