nextflow.enable.dsl=2

// workDir = '${params.rootDir}/workflows'
targetDir = "${params.rootDir}/target/nextflow"

include  { cellranger_mkfastq }       from  targetDir + "/demux/cellranger_mkfastq/main.nf"          params(params)
include  { cellranger_count }         from  targetDir + "/alignment/cellranger_count/main.nf"        params(params)


/* CellRanger - common workflow
 * 
 * consumed params:
 *   input:                         a path to input fastq files
 *   output                         a publish dir for the output
 * input format:                [ input, params ]
 *   value input:                   The input file directory
 *   value params:                  the params object
 * output format:               [ id, out_dir ]
 *   value id:                      same as input
 *   value out_dir:                 output directory
 * publishes:
 *   the output aggregated files
 */

// Nextflow entry point
workflow {
    
    main:
    // Auto-assigns ID
    if (!params.containsKey("id") || params.id == "") {
        id = "unk"
    } else {
        id = params.id
    }
    
    if (params.containsKey("input") && params.containsKey("sample_sheet") && params.containsKey("transcriptome")) {
        output_ch = Channel.value(
                [
                    id,
                    [
                        "input": file(params.input),
                        "sample_sheet": file(params.sample_sheet),
                    ],
                    params, 
                ]
            )
            | preprocessing_workflow
    } else {
        exit 1, "ERROR: Please provide --input, transcriptome, and --sample_sheet parameters"
    }

    emit: output_ch
}

workflow preprocessing_workflow {
    
    take:
    input_ch

    main:
    // Run count component
    output_ch = input_ch
        | cellranger_mkfastq
        | map { id, data, params -> [ id, [ input: file(data), transcriptome: file(params.transcriptome) ], params ]}
        | cellranger_count

    emit: output_ch
}