//
// Utility functions for the pipeline
//

class Utils {

    //
    // Function to generate software versions YAML
    //
    public static def softwareVersionsToYAML(ch_versions) {
        return ch_versions
                .unique()
                .map { it.text }
    }
}
