
// Prepare a manfiest
process prepare_manifest {
    conda params.tg_conda_env
    storeDir "${params.rn_image_dir}/${plate}"
    
    input:
        tuple val(plate), val(index_xml)
    output:
        path "manifest.tsv"
    script:
    """
    python $params.tg_core_dir/parse_xml.py \
    --input_file $index_xml \
    --output_path ./ \
    --to_manifest
    
    head -n 3 manifest.tsv > tmp
    mv tmp manifest.tsv
    
    """
    //manifest = params.tglow_image_dir + "/" + plate + "/" + "manifest.tsv"
}

// Fetches raw data from NFS and recodes into new OME file structure
// imaging queue
process fetch_raw {
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
    conda params.tg_conda_env
    publishDir "${params.rn_publish_dir}/basicpy/${plate}", mode: 'copy'

    input:
        tuple val(plate), val(img_channel)
    output:
        path "${plate}_ch${img_channel}"      
    script:
    """
    python $params.tg_core_dir/train_basicpy_2.py \
    --input $params.rn_image_dir \
    --output ./ \
    --output_prefix $plate \
    --plate $plate \
    --nimg 3 \
    --no_tune \
    --channel $img_channel
    """
}

// Cellpose
// gpu_normal queue
process cellpose {
    conda params.tg_conda_env
    scratch true
    storeDir "${params.rn_publish_dir}/masks/${plate}"
    
    input:
        tuple val(well), val(plate)
    script:
    """
    """
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
        val flatfields, optional: true
        val registration_matrices, optional: true
    output:
        path "*.tsv"
        
    script:
    
        // Outputs the cp files into ./images
        command = "\
        python stage_cellprofiler.py \
        --input $input_dir \
        --output ./images \
        "
        
        if (params.rn_max_project) {
            command = command + "--max_project "
        }
        
        // Run cell profiler
        command = command + "\n"
        command = command + "\
        cellprofiler \
        -c \
        -r \
        -p $params.cp_pipeline \
        -o ./ \
        -i ./images \
        --plugins-directory $params.cp_plugins"
        
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
        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml, tuple(row.channels.split(',')),  tuple(row.bp_channels.split(',')), row.cp_nucl_channel, row.cp_other_channel)}
            
        //manifest.view()

        if (params.bp_channels == null) {
            def csvFile = new File(params.rn_manifest)
            def csvData = csvFile.readLines()

            def plate_channel = []
            // Skip header
            // Substract one from the channel is python is 0 indexed
            for (int i = 1; i < csvData.size(); i++) {
                def curLine = csvData[i].split('\t')
                def curChannels = curLine[2].split(',')
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
        
        // Loop over previously generated manifests assuming stage has been run
        if (params.rn_manifest_well == null) {
            // code that reads paths available manfiests from previous stage
            manifests_in = Channel.fromPath("${params.rn_image_dir}/*/manifest.tsv")
        } else {
            manifests_in = Channel.from(params.rn_manifest_well)
        }

        //manifests_in.view()
        well_channel = manifests_in.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }

        // Re-order the well channel
        well_channel = well_channel.map{ row -> tuple(row.plate, row.well, row.row, row.col)}
        
        combined = well_channel.combine(manifest, by: 0)
        .map{ row -> tuple(plate: row[0], well: row[1], nucl_channel: row[7], other_channel: row[8])}
        .collect()
        .collect()

        combined.view()
        //combined.view()
        //combined.view()        
        

    // Cellpose - GPU
    // - plate
    // - well
    // - [OPTIONAL] estimated cell size
    // - [OPTIONAL] max project
    // + masks [publish?]

    // Register
    // - plate ref
    // - well
    // - second plates (1,2,3,4,5)
    // + registration matrices [publish]

    // Deconvelute - GPU
    // if params.decon:
        // - channel
        // - plate
        // - well
        // - [OPTIONAL] flatfield (perhaps not needed) 
        // img_channel = decon
    // else:
        // img_channel = raw


    // Cellprofiler
        // Stage
        // - plate
        // - well
        // - [OPTIONAL] max project
        // - [OPTIONAL] flatfield
        // - 

        // Cellprofiler
        // - image dir
        // - pipeline
        // + output files [export]
}

