## Get adjusted p-value and logFC matrices from SCE object
.getTestCols <- function(sce, adjp_cutoff = 0.0005, logfc_cutoff = 1) {
    rd <- rowData(sce)
    ## Get t-statistics and adjusted p-values from all comparisons
    tstat_cols <- setdiff(grep("\\.t$", colnames(rd), value = TRUE),
                          c(grep("highsaltcompl", colnames(rd), value = TRUE),
                            grep("atf1_150_tube|Atf1_150_tube", colnames(rd), value = TRUE)))
    adjp_cols <- setdiff(grep("adj.P.Val", colnames(rd), value = TRUE),
                         c(grep("highsaltcompl", colnames(rd), value = TRUE),
                           grep("atf1_150_tube|Atf1_150_tube", colnames(rd), value = TRUE)))
    pval_cols <- setdiff(grep("P.Value", colnames(rd), value = TRUE),
                         c(grep("highsaltcompl", colnames(rd), value = TRUE),
                           grep("atf1_150_tube|Atf1_150_tube", colnames(rd), value = TRUE)))
    logfc_cols <- setdiff(grep("logFC", colnames(rd), value = TRUE),
                          c(grep("highsaltcompl", colnames(rd), value = TRUE),
                            grep("se.logFC", colnames(rd), value = TRUE),
                            grep("atf1_150_tube|Atf1_150_tube", colnames(rd), value = TRUE)))

    tstat <- rd[, tstat_cols]
    adjp <- rd[, adjp_cols]
    pval <- rd[, pval_cols]
    logfc <- rd[, logfc_cols]

    stopifnot(all(sub("adj.P.Val", "t", colnames(adjp)) == colnames(tstat)))
    stopifnot(all(sub("adj.P.Val", "logFC", colnames(adjp)) == colnames(logfc)))
    stopifnot(all(sub("adj.P.Val", "P.Value", colnames(adjp)) == colnames(pval)))

    tstat <- as.matrix(tstat)
    adjp <- as.matrix(adjp)
    pval <- as.matrix(pval)
    logfc <- as.matrix(logfc)

    ## Define interactors based on the provided cutoffs
    interactor <- adjp <= adjp_cutoff &
        logfc > logfc_cutoff
    colnames(interactor) <- sub("\\.adj.P.Val", "", colnames(interactor))

    return(list(tstat = tstat, adjp = adjp, logfc = logfc,
                pval = pval, interactor = interactor))
}

