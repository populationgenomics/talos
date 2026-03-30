// not currently active, unsure if we want to pursue this
process ConvertSpliceVarDb {
    container params.container

    input:
        path(svdb)

    output:
        path("svdb.ht.tar.gz")

	script:
		"""
		set -euo pipefail

		python -m talos.convert_splice_var_db \
			--input ${svdb} \
			--output svdb.ht
		tar --no-xattrs -czf svdb.ht.tar.gz svdb.ht
		"""
}
