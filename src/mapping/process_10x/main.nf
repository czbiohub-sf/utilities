srcDir = params.rootDir + "/src"
targetDir = params.rootDir + "/module_openpipeline/target/nextflow"

include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { cellbender_remove_background } from targetDir + "/correction/cellbender_remove_background/main.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { add_id } from targetDir + "/metadata/add_id/main.nf"

include { processConfig; paramExists; readConfig; viashChannel; helpMessage; paramsToList } from srcDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from srcDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")
auto_config = readConfig("$projectDir/auto.vsh.yaml")

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
      auto: [ 
        publish: true
      ]
    )

    // rename .obs_names and add .obs["sample_id"]
    | pmap{ id, file ->
      def new_data = [
        input_id: id, 
        input: file, 
        obs_output: 'sample_id', 
        make_observation_keys_unique: true
      ]
      [ id, new_data ]
    }
    | add_id

    | pmap{ id, out, orig_data -> 
      [ id, out ]
    }

  emit: output_ch
}

/*
 * Integration test workflow.
 */
workflow test_wf {
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
}

/*
Workflow for automatically generating a `--param_list` YAML file.

Example. When providing this workflow with a directory with the structure
listed under "Input", the output of the pipeline will be similar to
that of the "Output".

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

import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.representer.Representer
import java.beans.IntrospectionException
import org.yaml.snakeyaml.introspector.Property

class NonMetaClassRepresenter extends Representer {
  protected Set<Property> getProperties( Class<? extends Object> type ) throws IntrospectionException {
    super.getProperties( type ).findAll { it.name != 'metaClass' }
  }
}

def writeParams(pl) {
  plStrings = pl.collect { data ->
    data.collectEntries{key, val -> 
      if (val instanceof List) {
        new_val = val.collect{it.toString()}
      } else {
        new_val = val.toString()
      }

      [key, new_val]
    }
  }
  def publishDir = 
    params.containsKey("publish_dir") ? params.publish_dir : 
    params.containsKey("publishDir") ? params.publishDir : 
    null
  def paramsFile = new File("${publishDir}/params.yaml")
  paramsFile.getParentFile().mkdirs()
  yaml = new Yaml(new NonMetaClassRepresenter())
  output = yaml.dump(plStrings)
  paramsFile.write(output)
}

workflow auto {
  helpMessage(auto_config)

  autoParams = paramsToList(params, auto_config)

  wfParams = autoParams.collectMany{ data ->
    def input_dir = file(data.input_dir)
    def fastq_glob = data.fastq_glob

    // look for 10x fastq 
    println("before glob: ${input_dir}/${fastq_glob}")
    def fastq_files = file("${input_dir}/${fastq_glob}")
    println("fastq_files: $fastq_files")

    // group by id
    // use glob to search for the sample id
    def sample_id_regex = fastq_glob.replace("**", "(.*)")
    def fastq_grouped = fastq_files.groupBy{ fastq_file ->
      fastq_file.toString()
        .replace(input_dir.toString() + "/", "") // strip parent directory
        .replaceAll(sample_id_regex, '$1') // extract sample id using regex
    }

    // create output list
    fastq_grouped.collect{ sample_id, input ->
      def output_raw = "${sample_id}_output_raw"
      def output_h5mu = "${sample_id}_output.h5mu"

      [
        id: sample_id,
        input: input,
        reference: data.reference,
        output_raw: output_raw,
        output_h5mu: output_h5mu
      ]
    }
  }
  // Log params file to output dir
  // todo: make sure they don't overwrite eachother
  writeParams(wfParams)

  Channel.fromList(wfParams)
    | map{tup -> [tup.id, tup]}
    | run_wf
}