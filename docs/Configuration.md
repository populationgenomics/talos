# Configuration File

The configuration file is a TOML file, containing headed sections for each principal step of Talos:

1. `GeneratePED` - no config section
2. `GeneratePanelData` - relevant to identifying a PanelApp instance to use, default panel to use, panels to add (this concerns both `DownloadPanelApp` and `UnifiedPanelAppParser`)
3. `RunHailFiltering` - relevant to both RunHailFiltering and RunHailFilteringSV, threshold values for various category tests, transcript consequences which are treated as `critical`, etc.
4. `ValidateMOI` - AF thresholds to apply during MOI tests, categories which are only retained in phenotype-matched genes, categories to remove from this analysis run
5. `CreateTalosHTML` - variables here are used in the CPG HTML implementation, mostly to connect internal and Seqr samples IDs, and to identify an appropriate seqr instance
6. `categories` - not Stage specific, but this contains a plaintext description of each category in the current codebase

This config file should be created in full prior to running Talos stages, a template is present at [example_config.toml](../src/talos/example_config.toml)
