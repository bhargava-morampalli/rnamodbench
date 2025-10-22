//
// This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
//

import org.everit.json.schema.Schema
import org.everit.json.schema.loader.SchemaLoader
import org.everit.json.schema.ValidationException
import org.json.JSONObject
import org.json.JSONTokener
import org.json.JSONArray
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

class NfcoreSchema {

    //
    // Function to loop over all parameters defined in schema and check
    // whether the given parameters adhere to the specifications
    //
    public static void validateParameters(params, json_schema, log) {
        // Simply return if no schema is provided
        if(json_schema == null) {
            return
        }
        def has_error = false
        //==================================================================//
        // Check for params not in the schema
        //
        params.each { param_name, param_value ->
            if(!json_schema.get('definitions').keySet().contains(param_name)) {
                if(!params.schema_ignore_params.toString().split(',').contains(param_name)) {
                    log.warn "Parameter '$param_name' is not recognised. Ignoring."
                }
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
