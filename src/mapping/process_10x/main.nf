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

/*
Recreate same folder structure:

Input:
/hpc/projects/data_lg/tabula_sapiens/renamed_fastqs/TSP1/fastqs/10X/
├── TSP1_Bladder_NA_10X_1_1
│   ├── TSP1_bladder_1_S13_L001_R1_001.fastq.gz
│   └── TSP1_bladder_1_S13_L001_R2_001.fastq.gz
└── foo
    └── TSP1_Pancreas_exocrine_10X_2_3
        ├── TSP1_exopancreas2_3_S12_L001_I1_001.fastq.gz
        └── TSP1_exopancreas2_3_S12_L001_R1_001.fastq.gz

Output:
/hpc/projects/data_lg/tabula_sapiens/realignment_gencode_v41/TSP1/mapping/10X/
├── TSP1_Bladder_NA_10X_1_1
│   ├── star_output
│   └── dataset.h5mu
└── foo
    └── TSP1_Pancreas_exocrine_10X_2_3
        ├── star_output
        └── dataset.h5mu
*/
workflow auto {
  helpMessage(auto_config)
  viashChannel(params, auto_config)
    | view{"original inputs: $it"}
    | flatMap{ id, data ->
      def input_dir = file(data.input_dir)
      // look for 10x fastq files
      def fastq_files = file("${input_dir}/**_S[0-9]*_L[0-9]*_R[12]_001.fastq.gz")

      // group by id
      def fastq_grouped = fastq_files.groupBy{ fastq_file ->
        fastq_file.toString().replace(input_dir.toString() + "/", "").replaceAll("_S[0-9]+_L[0-9]+_R[12]_001.fastq.gz", "")
      }

      // create output list
      fastq_grouped.collect{ sample_id, input ->
        def output_raw = "${sample_id}_output_raw"
        def output_h5mu = "${sample_id}_output.h5mu"

        def new_data = [
          input: input,
          reference: data.reference,
          output_raw: output_raw,
          output_h5mu: output_h5mu
        ]
        [sample_id, new_data]
      }
    }
    | view { id, data ->
      "$id:" +
        "  data.input.size(): ${data.input.size()}" +
        "  data.input[0]: ${data.input[0]}" +
        "  data.reference: ${data.reference}" +
        "  data.output_raw: ${data.output_raw}" +
        "  data.output_h5mu: ${data.output_h5mu}"
    }
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