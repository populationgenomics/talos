process GetLatestClinvArbitrationFile {
    publishDir params.large_files, mode: 'copy'

    container params.container

    input:
        val current_zenodo

    output:
        path "clinvarbitration_${current_zenodo}.tar.gz"

    shell:
    '''
    download_link=$(python -m talos.check_clinvarbitration_zenodo !{current_zenodo} --download)
    wget "${download_link}" -O clinvarbitration_!{current_zenodo}.tar.gz`
    '''
}
