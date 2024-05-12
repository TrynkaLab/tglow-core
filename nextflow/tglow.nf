
// Prepare a manfiest
process prepare_manifest {
    label 'tiny_img'
    conda params.tg_conda_env
    storeDir "${params.rn_image_dir}/${plate}"
    
    input:
        tuple val(plate), val(index_xml)
    output:
        path "manifest.tsv"
    script:
        cmd =
        """
        python $params.tg_core_dir/parse_xml.py \
        --input_file $index_xml \
        --output_path ./ \
        --to_manifest
        """
        
        // Only keep the first few lines 
        if (params.rn_testmode) {
            cmd += 
            """
            head -n 3 manifest.tsv > tmp
            mv tmp manifest.tsv
            """
        }
        cmd
    
    //manifest = params.tglow_image_dir + "/" + plate + "/" + "manifest.tsv"
}

// Fetches raw data from NFS and recodes into new OME file structure
// imaging queue
process fetch_raw {
    label 'tiny_img'
    conda params.tg_conda_env
    //storeDir "$params.rn_image_dir/$plate/$row/$col", mode: 'move'
    storeDir "${params.rn_image_dir}"

    input:
        tuple val(well), val(row), val(col), val(plate), val(index_xml)
    output:
        path "$plate/$row/$col"
    script:
        """
        python $params.tg_core_dir/convert_pe_raw.py \
        --input_file '$index_xml' \
        --output_path ./ \
        --well $well
        """         
}

// Basicpy
// normal queue
process basicpy {
    label 'himem'
    conda params.tg_conda_env
    publishDir "${params.rn_publish_dir}/basicpy/${plate}", mode: 'copy'

    input:
        tuple val(plate), val(img_channel)
    output:
        path "${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(plate), val(img_channel), emit: basicpy_out       
    script:
        cmd =
        """
        python $params.tg_core_dir/run_basicpy.py \
        --input $params.rn_image_dir \
        --output ./ \
        --output_prefix $plate \
        --plate $plate \
        --nimg 3 \
        --no_tune \
        --channel $img_channel \
        """
        
        if (params.rn_max_project) {
            cmd += "--max_project"
        }
            
        cmd
}

// Cellpose
process cellpose {
    label 'gpu_midmem'
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/masks/"
    //scratch true

    input:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel)
    output:
        path "${plate}/${row}/${col}/*_cell_mask*.tif", emit: cell_masks
        path "${plate}/${row}/${col}/*_nucl_mask*.tif", emit: nucl_masks, optional: true
        tuple val(plate), val(well), val(row), val(col), path("${plate}/${row}/${col}/*_cell_mask*.tif"), path("${plate}/${row}/${col}/*_nucl_mask*.tif"), emit: cellpose_out
    script:
        cmd =
        """
        python $params.tg_core_dir/run_cellpose.py \
        --input $params.rn_image_dir \
        --output ./ \
        --plate $plate \
        --well $well \
        --cell_channel $cell_channel \
        --gpu \
        --diameter $params.cp_cell_size \
        """
        
        if ($nucl_channel > 0) {
            cmd +=
            """
            --nucl_channel $nucl_channel  \
            --diameter_nucl $params.cp_nucl_size
            """
        }    
        
        if (params.rn_max_project) {
            cmd += "--no_3d"
        }
        cmd

}


// Run a cellprofiler run
// regular queue
process cellprofiler {
    scratch true
    conda params.cp_conda_env
    publishDir "$params.rn_publish_dir", mode: 'move'
    
    input:
        path input_dir
        val plate
        val well
        path "images/*" masks
        path flatfields, optional: true
        path registration_matrices, optional: true
    output:
        path "*.tsv"
        
    script:
    
        // Outputs the cp files into ./images
        command = 
        """
        
        # Stage files
        python stage_cellprofiler.py \
        --input $input_dir \
        --output ./images \
        """
        
        if (params.rn_max_project) {
            command += "--max_project"
        }
        
        // Run cell profiler
        command +=
        """
        # Run cellprofiler
        cellprofiler \
        -c \
        -r \
        -p $params.cp_pipeline \
        -o ./ \
        -i ./images \
        --plugins-directory $params.cp_plugins
        """
        
        return command
}


// Workflow to stage the data from NFS to lustre
workflow stage {
    main:
        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml)}
        
        //manifest.view()
        
        if (params.rn_manifest_well == null) {
            manifests_in = prepare_manifest(manifest)
        } else {
            manifests_in = Channel.from(params.rn_manifest_well)
        }

        //manifests_in.view()
        well_channel = manifests_in.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }
        
        // -> Channel.fromPath(manifest_path)
        //.splitCsv(header:false, sep:"\t")
        //.map( row -> tuple(val(row[0]), val(row[1]), val[row[2]]))
       // }
        
        //well_channel.view()
        
        fetch_raw(well_channel)
        
        // Read manifest
        //manifest = Channel.fromPath(params.manifest)
        //    .splitCsv(header:true)
        //    .map { row -> tuple(row.well, row.plate, row.index_xml)}
        
        // Fetch raw data
        //fetch_raw(manifest)       
}

// Main workflow
workflow run_pipeline {

    main:
    
        //------------------------------------------------------------
        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml, tuple(row.channels.split(',')),  tuple(row.bp_channels.split(',')), row.cp_nucl_channel, row.cp_cell_channel)}
            
        //manifest.view()

        //------------------------------------------------------------
        // Run basicpy
        
        // If there is no global overide on basicpy channels, get them from manifest
        if (params.bp_channels == null) {
            def csvFile = new File(params.rn_manifest)
            def csvData = csvFile.readLines()

            def plate_channel = []
            // Skip header
            // Substract one from the channel is python is 0 indexed
            for (int i = 1; i < csvData.size(); i++) {
                def curLine = csvData[i].split('\t')
                def curChannels = curLine[3].split(',')
                for (channel in curChannels) {
                    plate_channel << tuple(curLine[0], channel.toInteger()-1)
                }
            }
            basicpy_in = Channel.from(plate_channel)
        } else {
            basicpy_in = Channel.from(params.bp_channels)
        }
        
        //basicpy_channels = basicpy_channels.combine(manifest, by:0)
        //basicpy_in.view()
        basicpy_out = basicpy(basicpy_in)
        
        //------------------------------------------------------------
        // Loop over previously generated manifests assuming stage has been run
        // Here we start to run the seqeuntial processes, the first part runs independently
        if (params.rn_manifest_well == null) {
            // code that reads paths available manfiests from previous stage
            manifests_in = Channel.fromPath("${params.rn_image_dir}/*/manifest.tsv")
        } else {
            manifests_in = Channel.from(params.rn_manifest_well)
        }

        well_channel = manifests_in.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }

        // Re-order the well channel for later merging
        well_in = well_channel.map{ row -> tuple(row.plate, row.well, row.row, row.col)}
        
        //well_in.view()
        //------------------------------------------------------------
        // Cellpose
        
        // Append the things from the manifest needed for cellpose
        // Use -1 as the channels are 0 indexed
        cellpose_in = well_in.combine(manifest, by: 0)
        .flatMap{ row -> tuple(plate: row[0], well: row[1], row: row[2], col: row[3], nucl_channel: row[7]-1, cell_channel: row[8]-1)}

        //cellpose_in.view()
        
        cellpose_out = cellpose(cellpose_in)
        
        //------------------------------------------------------------
        // Deconvelute
        //deconvelute_out = deconvelute(well_in)
        
        //------------------------------------------------------------
        // Register
        // remap well in channel that combines plates based on a config item
        // so:
        // [plate 1, well, (ch1, ch2)]
        // [plate 2, well, (ch4, ch2)]
        // [plate N, well, (ch5, ch3)]

        // becomes:
        // [plate 1, well, ch1, ch2, [plate 2, plate N], [ch4, ch5], [ch2, ch3]]
        
        // Either output the same channel and do the actual merge later
        // just the registration matrices are saved, or save the actual registered
        // images
        
        //registration_out = register(well_in)
        
        
        //------------------------------------------------------------
        // Cellprofiler / get_features
        
        // Input channel:
        // plate, well, row, col, cellpose_nucl, cellpose_cell, registration_mat[optional], basicpy_dict[optional]
        
        // 

}

