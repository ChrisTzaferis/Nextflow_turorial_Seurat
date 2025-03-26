include { RUN_SEURAT } from './modules/local/runSeurat/main.nf'

workflow {
    ch_input = channel.fromPath(params.input_folder)
    ch_output_name = channel.value(params.output_name)
    ch_mito_pattern = channel.value(params.mt_pattern)

    // Capture the single output channel from RUN_SEURAT
    myOut = RUN_SEURAT(ch_input, ch_output_name, ch_mito_pattern)
    myOut.view()
}