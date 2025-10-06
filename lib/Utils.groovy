//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml
import java.nio.file.NoSuchFileException
import java.io.IOException
import nextflow.Nextflow

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

        /**
     * Resolve input paths with optional existence checking using idiomatic Nextflow patterns
     * 
     * @param log Nextflow logger
     * @param meta Sample metadata containing id and other identifiers  
     * @param base_dir Base directory path (can be local or S3)
     * @param relative_path Relative path from base_dir to target file
     * @param description Human-readable description for logging
     * @param optional Whether the file is optional (default: false)
     * @return Resolved path string or null for missing optional files
     */
    static def resolve_input_path(log, meta, base_dir, relative_path, description, boolean optional = false) {
        def resolved_path = Nextflow.file(base_dir).resolve(relative_path).toUriString()

        try {
            Nextflow.file(resolved_path, checkIfExists: true)
            return resolved_path
        } catch (NoSuchFileException e) {
            if (optional) {
                log.warn "Optional ${description} missing for sample ${meta.id} at ${resolved_path} - pipeline will continue without this file"
                return []
            } else {
                log.error "${description} not found for ${meta.id}: ${resolved_path}"
                Nextflow.exit(1)
            }
        }
    }

}
