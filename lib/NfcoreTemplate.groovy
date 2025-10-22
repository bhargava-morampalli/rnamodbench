//
// This file holds several functions used within the nf-core pipeline template
//

import org.yaml.snakeyaml.Yaml

class NfcoreTemplate {

    //
    // Check AWS Batch related parameters have been specified correctly
    //
    public static void awsBatch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            // Check params.awsqueue and params.awsregion have been set if running on AWSBatch
            assert (params.awsqueue && params.awsregion) : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            assert params.outdir.startsWith('s3:')       : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
        }
    }

    //
    // Print nf-core logo
    //
    public static String logo(workflow, monochrome_logs=true) {
        Map colors = logColours(monochrome_logs)
        String.format(
            """\n
            ${dashedLine(monochrome_logs)}
            ${colors.green}${workflow.manifest.name} v${workflow.manifest.version}${colors.reset}
            ${dashedLine(monochrome_logs)}
            """.stripIndent()
        )
    }

    //
    // Return dashed line
    //
    public static String dashedLine(monochrome_logs=true) {
        Map colors = logColours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }

    //
    // ANSII colours used for terminal logging
    //
    public static Map logColours(monochrome_logs=true) {
        Map colorcodes = [:]

        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        return colorcodes
    }

    //
    // Does what is says on the tin
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }

    //
    // Print a warning
    //
    public static void warn(log, message) {
        Map colors = logColours(true)
        log.warn("${colors.yellow}WARNING:${colors.reset} " + message)
    }

    //
    // ANSII Colours used for terminal logging
    //
    public static String logColours(Boolean monochrome_logs=true) {
        Map colorcodes = [:]

        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"

        return colorcodes
    }
}
