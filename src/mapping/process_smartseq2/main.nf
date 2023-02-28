nextflow.enable.dsl=2

srcDir = params.rootDir + "/src"
targetDir = params.rootDir + "/module_openpipeline/target/nextflow"

include { multi_star } from targetDir + "/mapping/multi_star/main.nf"
include { multi_star_to_h5mu } from targetDir + "/mapping/multi_star_to_h5mu/main.nf"

include { processConfig; paramExists; readConfig; viashChannel; helpMessage; paramsToList } from srcDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from srcDir + "/wf_utils/DataflowHelper.nf"

configDir = "${params.rootDir}/src/mapping/process_smartseq2"
config = readConfig("${configDir}/config.vsh.yaml")
auto_config = readConfig("${configDir}/auto.vsh.yaml")

/*
 * Main CLI workflow. See `nextflow run main.nf --help` documentation and usage.
 */
workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

/*
 * Main workflow.
 * 
 * This is used by the default workflow, integration workflow and the auto workflow.
 */
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
    | multi_star.run(
      auto: [publish: true]
    )
    | pmap{ id, output_dir, orig_data ->
      def new_data = [
        "input": output_dir,
        "output": orig_data.output_h5mu,
      ]
      [ id, new_data ]
    }
    | multi_star_to_h5mu.run(
      auto: [publish: true]
    )

  emit: output_ch
}

/*
 * Integration test workflow.
 */
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


/*
Workflow for automatically generating a `--param_list` YAML file.

Example. When providing this workflow with a directory with the structure
listed under "Input", the output of the pipeline will be similar to
that of the "Output".

Input:
/hpc/archives/AWS/buckets/tabula-sapiens/Pilot10/fastqs/smartseq2/
├── batch1
│   ├── TSP10_Fat_SCAT_SS2_B134180_B133833_Immune_N19_L002_R1.fastq.gz
│   └── TSP10_Fat_SCAT_SS2_B134180_B133833_Immune_N19_L002_R2.fastq.gz
├── TSP2_LungNeuron_proximal_SS2_B113453_B133092_Empty_G24_S288_R1_001.fastq.gz
└── TSP2_LungNeuron_proximal_SS2_B113453_B133092_Empty_G24_S288_R2_001.fastq.gz

Output:
/hpc/projects/data_lg/tabula_sapiens/realignment_gencode_v41/TSP1/mapping/10X/
├── batch1
│   └── TSP10_Fat_SCAT_SS2_B134180_B133833_Immune
│       ├── star_output
│       └── dataset.h5mu
└── TSP2_LungNeuron_proximal_SS2_B113453_B133092_Empty
    ├── star_output
    └── dataset.h5mu
*/

class NonMetaClassRepresenter extends org.yaml.snakeyaml.representer.Representer {
  protected Set<org.yaml.snakeyaml.introspector.Property> getProperties( Class<? extends Object> type ) throws java.beans.IntrospectionException {
    super.getProperties( type ).findAll { it.name != 'metaClass' }
  }
}

def getPublishDir() {
  def publishDir = 
    params.containsKey("publish_dir") ? params.publish_dir : 
    params.containsKey("publishDir") ? params.publishDir : 
    null
  return publishDir
}

def writeParams(param_list, params_file) {
  // convert file to strings
  param_list_strings = param_list.collect { data ->
    data.collectEntries{key, val -> 
      if (val instanceof List) {
        new_val = val.collect{it.toString()}
      } else {
        new_val = val.toString()
      }

      [key, new_val]
    }
  }

  // convert to yaml
  yaml = new org.yaml.snakeyaml.Yaml(new NonMetaClassRepresenter())
  output = yaml.dump(param_list_strings)

  // create parent directory
  params_file.getParent().mkdirs()

  // write to file
  params_file.write(output)
}

workflow auto {
  helpMessage(auto_config)

  // fetch params
  auto_params = paramsToList(params, auto_config)[0]

  // look for ss2 fastqs
  fastq_files = file("${auto_params.input_dir}/**.fastq.gz")
    .findAll{it.toString()
    .matches(auto_params.fastq_regex)}

  // create output list

  input_r1 = fastq_files.findAll{it.toString().matches(".*_R1[_.].*")}.sort()
  input_r2 = fastq_files.findAll{it.toString().matches(".*_R2[_.].*")}.sort()
  input_id = input_r1.collect{ fastq_file ->
    fastq_file.toString()
      .replace("${auto_params.input_dir}/", "") // remove root directory
      .replaceAll(auto_params.fastq_regex, auto_params.cell_id_replacement) // extract sample id using regex
  }

  assert input_r1.size() == input_r2.size(): "The number of R1 files (${input_r1.size()}) should equal the number of R2 files (${input_r2.size()})."

  param_list = [[
    id: auto_params.id,
    input_id: input_id,
    input_r1: input_r1,
    input_r2: input_r2,
    reference_index: auto_params.reference_index,
    reference_gtf: auto_params.reference_gtf,
    output_raw: "${auto_params.id}_raw",
    output_h5mu: "${auto_params.id}.h5mu"
  ]]

  // Log params file to output dir
  param_file = file("${getPublishDir()}/param_list.yaml")
  writeParams(param_list, param_file)

  // run pipeline
  Channel.fromList(param_list)
    | map{tup -> [tup.id, tup]}
    | run_wf
}