// record the annotated shards produced for a cohort, plus the settings that produced them.
// subsequent runs consume this manifest instead of globbing the annotated directory - it both
// gates re-annotation (manifest exists == cohort is done) and names exactly the shards to reuse,
// so shards left behind by earlier runs (e.g. a different vcf_split_n) can never leak in.
// native (exec) process - runs on the head node, no container required
process WriteShardManifest {
    publishDir { annotated_dir }, mode: 'copy'

    input:
        tuple val(cohort), val(shard_names), val(source), val(type), val(annotated_dir)

    output:
        tuple val(cohort), path("${cohort}_manifest.json")

    exec:
        def manifest = [
            cohort: cohort,
            source: source.toString(),
            type: type,
            vcf_split_n: params.vcf_split_n,
            n_shards: shard_names.size(),
            shards: shard_names.sort(false),
        ]
        task.workDir.resolve("${cohort}_manifest.json").text =
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(manifest))
}
