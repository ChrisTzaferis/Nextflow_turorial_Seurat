process RUN_SEURAT {
    conda "${moduleDir}/r_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"
        
    input:
    val input_folder_path
    val output_name
    val reg_exp_mito

    output:
    path "${output_name}.RDS"

    script:
    """
    Rscript ${baseDir}/bin/seuratAnalysis.R "${input_folder_path}" "${output_name}" "${reg_exp_mito}"
    """
}
