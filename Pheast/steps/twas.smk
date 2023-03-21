rule twas_compute_weights:
    """Use FUSION to compute TWAS weights from expression and genotypes."""
    input:
    output:
    shell:
        """
        Rscript FUSION.compute_weights.R \
            --bfile $INP \
            --tmp $TMP \
            --out $OUT \
            --models top1,blup,bslmm,lasso,enet
        """
