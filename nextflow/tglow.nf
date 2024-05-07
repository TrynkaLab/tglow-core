

// Prepare a manfiest
process prepare_manifest {
    conda params.conda_env
    publishDir "$params.tglow_image_dir/$plate", mode: 'copy'
    
    input:
        tuple val(plate), val(index_xml)
    output:
        path "manifest.tsv"
    script:
    """
    python $params.tglow_core_dir/parse_xml.py \
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
    conda params.conda_env
    publishDir "$params.tglow_image_dir/$plate/$row/$col", mode: 'move'

    input:
        tuple val(well), val(row), val(col), val(plate), val(index_xml)
    output:
        path "$plate/$row/$col"
    script:
    """
    python $params.tglow_core_dir/convert_pe_raw.py \
    --input_file '$index_xml' \
    --output_path ./ \
    --well $well
    """         
}


workflow stage {
    main:
        // Read manifest
        manifest = Channel.fromPath(params.manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml)}
        
        //manifest.view()
        
        if (params.plate_manifest == null) {
            plate_channel = prepare_manifest(manifest)
        } else {
            plate_channel = Channel.from(params.plate_manifest)
        }

        //plate_channel.view()
        well_channel = plate_channel.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }
        
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