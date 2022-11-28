nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/modules_openpipeline/target/nextflow"

// include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { join_csv } from targetDir + "/metadata/join_csv/main.nf"
include { join_uns_to_obs } from targetDir + "/metadata/join_uns_to_obs/main.nf"
include { concat } from targetDir + "/dataflow/concat/main.nf"
include { from_h5mu_to_h5ad } from targetDir + "/convert/from_h5mu_to_h5ad/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataFlowHelper.nf"

config = readConfig("$workflowDir/concat/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}


workflow run_wf {
  take: input_ch

  main:
    output_ch = input_ch
      // rename .obs_names and add .obs["sample_id"]
      | pmap{ id, data ->
        def new_data = [
          input_id: id, 
          input: data.input, 
          obs_output: 'sample_id', 
          make_observation_keys_unique: true
        ]
        [ id, new_data, data ]
      }
      | add_id

      // join metadata csv to .obs
      | pmap{ id, file, orig_data ->
        def new_data = [ 
          input: file, 
          input_csv: orig_data.input_metadata,
          obs_key: 'sample_id'
        ]
        [ id, new_data, orig_data ]
      }
      | join_csv

      // temporary step to fix .obs
      | pmap{ id, file, orig_data -> 
        def new_data = [ 
          input: file,
          uns_key: "metrics_cellranger"
        ]
        [ id, new_data, orig_data ]
      }
      | join_uns_to_obs

      // combine into one channel event
      | toSortedList{ a, b -> b[0] <=> a[0] }
      | map { tups -> 
        def new_data = [ 
          input_id: tups.collect{it[0]}, 
          input: tups.collect{it[1]}
        ]
        [ "combined", new_data, tups[0][2] ]
      }

      // concatenate into one h5mu
      | pmap{ id, data, other ->
        [ id, data + [ output: other.output ], other]
      }
      | concat.run(
        auto: [ publish: true ]
      )

      // convert to h5ad
      | pmap{ id, file, other ->
        def output_h5ad = other.output.replaceAll('.h5mu$', ".h5ad")
        [ id, [ input: file, output: output_h5ad ], other]
      }
      | from_h5mu_to_h5ad.run(
        auto: [ publish: true ]
      )

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
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
        input_type: "10xh5",
        modality: null
      ],
      [
        id: "10xmtx_test",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix",
        input_type: "10xmtx",
        modality: null,
        output: "\$id.h5mu"
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