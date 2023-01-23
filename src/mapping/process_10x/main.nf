nextflow.enable.dsl=2

srcDir = params.rootDir + "/src"
targetDir = params.rootDir + "/module_openpipeline/target/nextflow"

include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { cellbender_remove_background } from targetDir + "/correction/cellbender_remove_background/main.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"

include { paramExists; readConfig; viashChannel; helpMessage } from srcDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from srcDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")
auto_config = readConfig("$projectDir/auto.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take: input_ch

  main:
    output_ch = input_ch

    // map data to reference
    | pmap{ id, data ->
      new_data = data + [ output: data.output_raw ]
      [ id, new_data, data ]
    }
    | cellranger_count.run(auto: [ publish: true ])

    // split output dir into map
    | cellranger_count_split

    // convert to h5mu
    | pmap{ id, data, orig_data -> 
      new_data = [ 
        input: data.raw_h5,
        input_metrics_summary: data.metrics_summary
      ]
      [ id, new_data, orig_data ]
    }
    | from_10xh5_to_h5mu

    // run cellbender
    | cellbender_remove_background.run(
      args: [
        expected_cells_from_qc: true,
        layer_output: "cellbender"
      ]
    )

    // filter counts
    | pmap{ id, file, orig_data -> 
      new_data = [ input: file, output: orig_data.output_h5mu ]
      [ id, new_data, orig_data ]
    }
    | filter_with_counts.run(
      args: [
        layer: "cellbender",
        min_genes: 100, 
        min_counts: 1000, 
        do_subset: true
      ],
      auto: [ publish: true ]
    )

    | pmap{ id, out, orig_data -> 
      [ id, out ]
    }

  emit: output_ch
}

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  testParams = [
    param_list: [
      [
        id: "pbmc_1k_v3",
        input: params.resources_test + "/pbmc_1k_v3/fastqs",
        reference: params.resources_test + "/reference_gencodev41_chr1/reference_cellranger.tar.gz"
      ]
    ]
  ]
    
  output_ch =
    viashChannel(testParams, config)
    | view { "Input: $it" }    
    | run_wf
    // | view { output ->
    //     assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    //     assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
    //     "Output: $output"
    //   }
    //   | toSortedList()
    //   | map { output_list ->
    //     assert output_list.size() == 3 : "output channel should contain three events"
    //   }
}