import Utils

// Test workflow for Utils.resolveInputPath function
workflow TEST_RESOLVE_INPUT_PATH {
    take:
        test_data // Channel of test data: [meta, base_dir, relative_path, description, optional]

    main:
        ch_results = test_data.map { meta, base_dir, relative_path, description, optional ->
            def result = Utils.resolveInputPath(meta, base_dir, relative_path, description, optional, log)
            return [meta, result]
        }

    emit:
        results = ch_results // Channel: [meta, resolved_path_or_null]
}
