nextflow.enable.dsl=2

srcDir = params.rootDir + "/src"
targetDir = params.rootDir + "/module_openpipeline/target/nextflow"

include { star_align } from targetDir + "/mapping/star_align/main.nf"
include { samtools_sort } from targetDir + "/mapping/samtools_sort/main.nf"
include { htseq_count } from targetDir + "/mapping/htseq_count/main.nf"
include { htseq_count_to_h5mu } from targetDir + "/mapping/htseq_count_to_h5mu/main.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"

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
        "input": data.input,
        "reference": data.reference_index,
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
        "outReadsUnmapped": "Fastx"
        // "genomeLoad": "LoadAndKeep",
      ]
      [ id, new_data, data ]
    }
    | star_align

    // sort aligned reads
    | pmap { id, dir, passthrough ->
      def new_data = [
        "input": workflow.stubRun ? dir : file(dir + "/Aligned.out.bam"),
        "output_bam": '$id.$key.output.bam',
        "output_bai": '$id.$key.output.bam.bai',
        "output_format": "bam"
      ]
      [id, new_data, passthrough]
    }
    | samtools_sort

    // create count matrix
    // args interpreted from: https://github.com/czbiohub/utilities/blob/47371ab0465d85ab8a2a4dfa092581385064ef62/src/utilities/alignment/run_star_and_htseq.py#L286
    | pmap { id, data, passthrough ->
      def new_data = [
        "input": data.output_bam,
        "reference": passthrough.reference_gtf,
        "order": "name",
        "stranded": "no",
        "output_sam_format": "bam",
        "mode": "intersection-nonempty"
      ]
      [id, new_data, passthrough]
    }
    | htseq_count
    
    // group events by sample id and convert to h5mu
    // todo: rename sample_id to group_id
    | map{ tup ->
      def id = tup[0]
      def data = tup[1]
      def passthrough = tup[2]
      def others = tup.drop(3)
      [ passthrough.sample_id, [input_id: id, input_counts: data.output, passthrough: passthrough, others: others]]
    }
    | groupTuple()
    | map{ id, tups ->
      def input_id = tups.collect{it.input_id}
      def input_counts = tups.collect{it.input_counts}
      def passthrough = tups[0].passthrough
      def others = tups[0].others
      def new_data = [
        "input_id": input_id,
        "input_counts": input_counts,
        "reference": passthrough.reference_gtf
      ]
      [ id, new_data, passthrough ] + others
    }
    | htseq_count_to_h5mu.run(
      auto: [ publish: true ]
    )

    // remove passthrough
    | pmap{ id, out, passthrough -> 
      [ id, out ]
    }

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