//
// This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
//

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

class NfcoreSchema {

    //
    // Function to loop over all parameters defined in schema and check
    // whether the given parameters adhere to the specifications
    //
    public static void validateParameters(workflow, params, log) {
        // If no schema is present, we can't validate
        def schema_path = "${workflow.projectDir}/nextflow_schema.json"
        if (!file(schema_path).exists()) {
            log.warn "Could not find pipeline schema file: ${schema_path}"
            return
        }

        // Parse the schema
        def schema = new JsonSlurper().parseText(file(schema_path).text)
        def schema_params = []
        schema.definitions.each { key, val ->
            if(val.properties) {
                val.properties.each { pkey, pval ->
                    schema_params.add(pkey)
                }
            }
        }

        // Check for parameters that are not in the schema
        def ignored_params = ['help', 'version', 'outdir', 'input', 'validate_params', 'monochrome_logs', 'schema_ignore_params']
        params.each { param_name, param_value ->
            if(!schema_params.contains(param_name) && !ignored_params.contains(param_name)) {
                log.warn "Parameter `${param_name}` is not recognised by the pipeline."
            }
        }
    }

    //
    // Function to get parameter information from the schema
    //
    public static Map paramsHelp(workflow, params, command) {
        return [:]
    }

    //
    // Function to print parameter summary log
    //
    public static Map paramsSummaryLog(workflow, params) {
        def summary_params = [:]
        summary_params['Core Nextflow options'] = ''
        summary_params['Input/output options'] = ''
        if (params.input) {
            summary_params['Input/output options'] += "  input: ${params.input}\n"
        }
        if (params.outdir) {
            summary_params['Input/output options'] += "  outdir: ${params.outdir}\n"
        }
        summary_params['Reference genome options'] = ''
        if (params.ref_16s) {
            summary_params['Reference genome options'] += "  ref_16s: ${params.ref_16s}\n"
        }
        if (params.ref_23s) {
            summary_params['Reference genome options'] += "  ref_23s: ${params.ref_23s}\n"
        }
        summary_params['Max job request options'] = ''
        summary_params['Max job request options'] += "  max_cpus: ${params.max_cpus}\n"
        summary_params['Max job request options'] += "  max_memory: ${params.max_memory}\n"
        summary_params['Max job request options'] += "  max_time: ${params.max_time}\n"

        return summary_params
    }

    //
    // Beautify parameters for summary and return as string
    //
    public static String paramsSummaryMap(workflow, params) {
        def summary_section = ''
        return summary_section
    }
}
