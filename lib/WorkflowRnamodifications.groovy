//
// This file holds several functions specific to the workflow/rnamodifications.nf in the nf-core/rnamodifications pipeline
//

class WorkflowRnamodifications {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (!params.input) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }

        if (!params.ref_16s) {
            log.error "Please provide a 16S rRNA reference file e.g. '--ref_16s 16s.fasta'"
            System.exit(1)
        }

        if (!params.ref_23s) {
            log.error "Please provide a 23S rRNA reference file e.g. '--ref_23s 23s.fasta'"
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }
}
