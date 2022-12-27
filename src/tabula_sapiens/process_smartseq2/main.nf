nextflow.enable.dsl=2

srcDir = params.rootDir + "/src"
targetDir = params.rootDir + "/module_openpipeline/target/nextflow"

include { multi_star } from targetDir + "/mapping/multi_star/main.nf"
include { multi_star_to_h5mu } from targetDir + "/mapping/multi_star_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from srcDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from srcDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take: input_ch

  main:
    output_ch = input_ch

    // map data to reference using STAR
    // args interpreted from https://github.com/czbiohub/utilities/blob/47371ab0465d85ab8a2a4dfa092581385064ef62/src/utilities/alignment/run_star_and_htseq.py#L29
    | pmap{ id, data ->
      def new_data = [
        "input_id": data.input_id,
        "input_r1": data.input_r1,
        "input_r2": data.input_r2,
        "reference_index": data.reference_index,
        "reference_gtf": data.reference_gtf,
        "output": data.output_raw,
        // STAR parameters
        "outFilterType": "BySJout",
        "outFilterMultimapNmax": 20,
        "alignSJoverhangMin": 8,
        "alignSJDBoverhangMin": 1,
        "outFilterMismatchNmax": 999,
        "outFilterMismatchNoverLmax": 0.04,
        "alignIntronMin": 20,
        "alignIntronMax": 1000000,
        "alignMatesGapMax": 1000000,
        "outSAMstrandField": "intronMotif",
        "outSAMtype": [ "BAM", "Unsorted" ],
        "outSAMattributes": [ "NH", "HI", "NM", "MD" ],
        "outReadsUnmapped": "Fastx",
        // HTSeq parameters
        "stranded": "no",
        "mode": "intersection-nonempty"
      ]
      [ id, new_data, data ]
    }
    | multi_star
    | pmap{ id, output_dir, orig_data ->
      def new_data = [
        "input": output_dir,
        "output": orig_data.output_h5mu,
      ]
      [ id, new_data ]
    }
    | multi_star_to_h5mu

  emit: output_ch
}

// workflow test_wf {
//   helpMessage(config)

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"

//   testParams = [
//     param_list: [
//       [
//         id: "pbmc_1k_v3",
//         input: params.resources_test + "/pbmc_1k_v3/fastqs",
//         reference: params.resources_test + "/reference_gencodev41_chr1/reference_cellranger.tar.gz"
//       ]
//     ]
//   ]
    
//   output_ch =
//     viashChannel(testParams, config)
//     | view { "Input: $it" }    
//     | run_wf
//     // | view { output ->
//     //     assert output.size() == 2 : "outputs should contain two elements; [id, file]"
//     //     assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
//     //     "Output: $output"
//     //   }
//     //   | toSortedList()
//     //   | map { output_list ->
//     //     assert output_list.size() == 3 : "output channel should contain three events"
//     //   }
// }