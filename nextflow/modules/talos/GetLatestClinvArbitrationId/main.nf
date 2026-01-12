process GetLatestClinvArbitrationId {
    container params.container

    input:
        val current_zenodo

    output:
        stdout

    """
    python -m talos.check_clinvarbitration_zenodo ${current_zenodo}
    """
}
