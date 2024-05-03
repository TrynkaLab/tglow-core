
// Fetches raw data from NFS and recodes into new OME file structure
// imaging queue
process fetch_raw {
    conda params.conda_env
    
    input:
        tuple val(well), val(plate), path(input_xml)
    script:
    """
    python $params.tglow_core_dir/convert_pe_raw.py \
    --input_file $input_xml \
    --output_path $params.tglow_image_dir \
    --well $well
    """         
}


workflow {
    
    // Read manifest
    manifest = Channel.fromPath(params.manifest)
        .splitCsv(header:true)
        .map { row -> tuple(row.well, row.plate, row.index_xml)}
    
    // Fetch raw data
    fetch_raw(manifest)    
        
}