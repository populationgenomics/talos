process CreateRoiFromGff3 {
    container params.container

    // this file is not downloaded at runtime - that presents issues for execution environments
    // disconnected from the internet
    input:
		path gff

    // generates a BED file one-gene-per-line, for per-gene annotations
    // generates an BED with overlapping intervals merged for filtering
    // sorts both
    output:
        path "GRCh38.bed", emit: bed
        path "GRCh38_merged.bed", emit: merged_bed

    script:
        """
        python -m talos.annotation_scripts.create_roi_from_gff3 \
            --gff3 ${gff} \
            --unmerged_output unsorted_GRCh38.bed \
            --merged_output unsorted_merged_GRCh38.bed

        sort -k1,1V -k2,2n unsorted_GRCh38.bed > GRCh38.bed
        sort -k1,1V -k2,2n unsorted_merged_GRCh38.bed > GRCh38_merged.bed
        """
}
