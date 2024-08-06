filterForHeatmap <- function(tstatmat, adjpmat, logfcmat, adjpThreshold,
                             logfcThreshold, maxNbrTargets) {
    ## Filter - require significant in at least one comparison
    tstatfilt <- tstatmat
    tstatfilt[adjpmat >= adjpThreshold] <- NA
    tstatfilt[logfcmat <= logfcThreshold] <- NA

    ## Some experiments pull down a large number of proteins
    ## -> remove proteins that are only being pulled down in these experiments
    keep <- which(colSums(!is.na(tstatfilt)) <= maxNbrTargets)
    tstatfilt <- tstatfilt[rowSums(!is.na(tstatfilt[, keep])) > 0, ]
    tstatfilt[is.na(tstatfilt)] <- 0
    tstatfilt <- tstatfilt[rowSums(tstatfilt) > 0, ]

    ## Subset the t-statistic matrix to only the retained proteins
    stopifnot(colnames(tstatfilt) == colnames(tstatmat))
    tstatmat <- tstatmat[rownames(tstatfilt), colnames(tstatfilt)]
    colnames(tstatmat) <- .getSimplifiedComparison(colnames(tstatmat))

    list(tstats = tstatmat, tstatsfilt = tstatfilt)
}

makeComplexAnnotation <- function(tstatmat, complexes, idmap,
                                  show_legend = TRUE) {
    rowannot <- data.frame(Complex = rep(NA, nrow(tstatmat)),
                           row.names = rownames(tstatmat))
    rowcols <- list(Complex = c())
    for (cmplx in names(complexes)) {
        rowannot$Complex[.getPomBaseIdFromProtein(rownames(rowannot),
                                                  idmap = idmap) %in%
                             complexes[[cmplx]]$pombase_ids] <- cmplx
        rowcols$Complex[cmplx] <- complexes[[cmplx]]$color
    }
    rowAnnotation(df = rowannot,
                  col = rowcols,
                  annotation_name_gp = gpar(fontsize = fl),
                  simple_anno_size_adjust = TRUE,
                  na_col = bg_color,
                  show_annotation_name = FALSE,
                  show_legend = show_legend,
                  annotation_legend_param = list(title_gp = gpar(fontsize = fl),
                                                 legend_direction = "horizontal"))
}

makeComplexAndDBDAnnotation <- function(tstatmat, complexes, dbd, idmap) {
    rowannot <- data.frame(Complex = rep(NA, nrow(tstatmat)),
                           row.names = rownames(tstatmat))
    rowcols <- list(Complex = c())
    for (cmplx in names(complexes)) {
        rowannot$Complex[.getPomBaseIdFromProtein(rownames(rowannot),
                                                  idmap = idmap) %in%
                             complexes[[cmplx]]$pombase_ids] <- cmplx
        rowcols$Complex[cmplx] <- complexes[[cmplx]]$color
    }
    rowannot$`DBD fam.` <- dbd$DBD_class_name[match(
        .getPomBaseIdFromProtein(rownames(tstatmat), idmap = idmap),
        dbd$PomBaseID)]
    setfam <- c("KilA-N", "Zn(II)2Cys6", "bZIP", "bHLH")
    unsetfam <- setdiff(unique(dbd$DBD_class_name), setfam)
    rowannot$`DBD fam.`[rowannot$`DBD fam.` %in% unsetfam] <- "Other (combined)"
    rowcols <- c(rowcols, list(`DBD fam.` = c(`KilA-N` = main_colors[5],
                                              `Zn(II)2Cys6` = main_colors[4],
                                              bZIP = main_colors[3],
                                              bHLH = main_colors[2],
                                              `Other (combined)` = bg_color)))
    rowAnnotation(df = rowannot,
                  col = rowcols,
                  annotation_name_gp = gpar(fontsize = fl),
                  simple_anno_size_adjust = TRUE,
                  na_col = bg_color,
                  show_annotation_name = TRUE,
                  show_legend = FALSE,
                  annotation_legend_param = list(title_gp = gpar(fontsize = fl)))
}

makeTFIntAnnotation <- function(baitdata) {
    HeatmapAnnotation(
        which = "column",
        `Number of\ninteractors` =
            anno_barplot(x = rowSums(as.matrix(baitdata$nbrIntMat[colnames(baitdata$mat),
                                                                  c("nbrTFs_samefam",
                                                                    "nbrTFs_difffam")])),
                         gp = gpar(fill = na_color,
                                   col = na_color),
                         border = FALSE, bar_width = 1.0,
                         axis_param = list(gp = gpar(fontsize = fs)),
                         width = unit(10, "mm")),
        annotation_name_side = "left",
        annotation_name_rot = 0,
        annotation_name_gp = gpar(fontsize = fl))
}

makeColLabels <- function(tstatmat, colLabel, fontsize) {
    columnAnnotation(
        TFnames = anno_mark(which = "column", side = "bottom",
                            at = match(colLabel, colnames(tstatmat)),
                            labels = .getProteinNameFromSimplifiedComparison(colLabel),
                            labels_gp = gpar(fontsize = fontsize)))
}

makeRowLabels <- function(tstatmat, rowLabel, fontsize) {
    rowAnnotation(
        TFnames = anno_mark(which = "row", side = "left",
                            at = match(rowLabel, rownames(tstatmat)),
                            labels = rowLabel,
                            labels_gp = gpar(fontsize = fontsize)))
}

makeHeatmapCol <- function(stringency = "high") {
    if (stringency == "high") {
        circlize::colorRamp2(
            breaks = seq(3.0, 13.0, length.out = 64),
            colors = colorRampPalette(binary_heatmap_colors[c("FALSE", "TRUE")])(64))
    } else {
        circlize::colorRamp2(
            breaks = seq(1.5, 13.0, length.out = 64),
            colors = colorRampPalette(binary_heatmap_colors[c("FALSE", "TRUE")])(64))
    }
}

reorderColsByRows <- function(tstat) {
    ## Row names should be protein names, column names simplified comparison names
    ord <- match(rownames(tstat), .getProteinNameFromSimplifiedComparison(colnames(tstat)))
    ord <- ord[!is.na(ord)]
    ord <- c(ord, setdiff(seq_len(ncol(tstat)), ord))
    tstat[, ord]
}

makeHeatmapData <- function(sce, adjpthr, log2fcthr, conc) {
    tc <- .getTestCols(sce, adjp_cutoff = adjpthr, logfc_cutoff = log2fcthr)

    tstats <- filterForHeatmap(tstatmat = tc$tstat, adjpmat = tc$adjp,
                               logfcmat = tc$logfc, logfcThreshold = log2fcthr,
                               adjpThreshold = adjpthr, maxNbrTargets = Inf)
    tstat <- tstats$tstats

    ## tstats$tstatsfile: protein x experiment
    ## Number of TFs that are pulled down by each bait
    all_tfs <- setdiff(.getProteinNameFromComparison(colnames(tstats$tstatsfilt)),
                       c("Untagged", "untagged", NA))
    nbrInt <- colSums(tstats$tstatsfilt[rownames(tstats$tstatsfilt) %in% all_tfs, ] > 0,
                      na.rm = TRUE)
    names(nbrInt) <- .getSimplifiedComparison(names(nbrInt))

    ## Split by family - return a matrix with several columns:
    ## #pulled down TFs (including itself), #pulled down TFs (excluding itself),
    ## #pulled down TFs from the same family (excluding itself),
    ## #pulled down TFs from another family
    tmptf <- tstats$tstatsfilt[match(.getProteinNameFromComparison(colnames(tstats$tstatsfilt)),
                                     rownames(tstats$tstatsfilt)), ]
    rownames(tmptf) <- .getProteinNameFromComparison(colnames(tstats$tstatsfilt))
    dbdtf <- data.frame(PomBaseID = .getPomBaseIdFromComparison(colnames(tstats$tstatsfilt),
                                                                idmap = idmap),
                        bait = .getProteinNameFromComparison(colnames(tstats$tstatsfilt)),
                        DBD_class_name = NA)
    for (i in seq_len(nrow(dbdtf))) {
        lookfor <- dbdtf$PomBaseID[i]
        if (!(lookfor %in% c("untagged", "Untagged"))) {
            dbdtf$DBD_class_name[i] <- dbd$DBD_class_name[match(lookfor, dbd$PomBaseID)]
        } else {
            dbdtf$DBD_class_name[i] <- "no_class"
        }
    }
    stopifnot(dbdtf$bait == .getProteinNameFromComparison(colnames(tmptf)))
    nbrIntMat <- do.call(dplyr::bind_rows, lapply(seq_along(colnames(tmptf)), function(i) {
        expr <- colnames(tmptf)[i]
        fam <- dbdtf$DBD_class_name[i]
        data.frame(experiment = expr,
                   bait = dbdtf$bait[i],
                   dbdfam = fam,
                   TFs = paste(rownames(tmptf)[-i][which(tmptf[-i, expr] > 0)],
                               collapse = ","),
                   nbrTFs = sum(tmptf[, expr] > 0, na.rm = TRUE),
                   nbrTFs_wobait = sum(tmptf[-i, expr] > 0, na.rm = TRUE),
                   nbrTFs_samefam = sum(tmptf[-i, expr] > 0 &
                                            dbdtf$DBD_class_name[-i] == fam &
                                            dbdtf$DBD_class_name[-i] != "Other",
                                        na.rm = TRUE),
                   nbrTFs_difffam = sum(tmptf[-i, expr] > 0 &
                                            (dbdtf$DBD_class_name[-i] != fam |
                                                 dbdtf$DBD_class_name[-i] == "Other"),
                                        na.rm = TRUE)
        )
    }))


    ## Split experiments into three groups
    stopifnot(all(colnames(tstat) == .getSimplifiedComparison(colnames(tstats$tstatsfilt))))
    if (conc == 150) {
        tstat_no_bait <-
            tstat[, colSums(tstats$tstatsfilt) == 0 |
                      (.getOrigBaitNameFromComparison(colnames(tstats$tstatsfilt), idmap = idmap) %in%
                           baitclass$Gene_name[baitclass$class == "nobait"])]
        tstat_150only <-
            tstat[, colSums(tstats$tstatsfilt) > 0 &
                      (.getOrigBaitNameFromComparison(colnames(tstats$tstatsfilt), idmap = idmap) %in%
                           baitclass$Gene_name[baitclass$class == "150only"])]
        tstat_150and500 <-
            tstat[, colSums(tstats$tstatsfilt) > 0 &
                      (.getOrigBaitNameFromComparison(colnames(tstats$tstatsfilt), idmap = idmap) %in%
                           baitclass$Gene_name[baitclass$class == "150and500"])]

        proteins_150only <- .getProteinNameFromSimplifiedComparison(colnames(tstat_150only))
        stopifnot(all(proteins_150only %in% rownames(tstat)))
        tstat_150and500 <- tstat_150and500[!(rownames(tstat_150and500) %in% proteins_150only), ]

        ## Cluster proteins based on the correlation between their t-stat vectors
        dst <- sqrt(2 * (1 - cor(t(tstat_150and500), use = "pairwise.complete")))
        dst[is.na(dst)] <- sqrt(2 * (1 - (-1)))
        clst <- hclust(as.dist(dst))

        ## Reorder proteins by the clustering
        protorder <- c(rownames(tstat_150and500)[clst$order], proteins_150only)
        tstat <- tstat[protorder, ]

        ## Reorder experiments to keep a "diagonal" in the heatmap
        tstat <- reorderColsByRows(tstat)

        tmpbaits <- .getOrigBaitNameFromSimplifiedComparison(colnames(tstat),
                                                             idmap = idmap)
        colsplit <- ifelse(
            tmpbaits %in% baitclass$Gene_name[baitclass$class == "nobait"],
            "No\nbait",
            ifelse(
                tmpbaits %in% baitclass$Gene_name[baitclass$class == "150only"],
                "150 mM\nNaCl IP-MS",
                ifelse(
                    tmpbaits %in% baitclass$Gene_name[baitclass$class == "150and500"],
                    "150 and 500 mM\nNaCl IP-MS", "NA"
                )
            )
        )
    } else if (conc == 500) {
        tstat_onlybait <-
            tstat[, colSums(tstats$tstatsfilt > 0) == 1 |
                      (.getOrigBaitNameFromComparison(colnames(tstats$tstatsfilt),
                                                      idmap = idmap) %in%
                           baitclass$Gene_name[baitclass$class500 == "Lost in 500 mM"])]
        tstat_several <-
            tstat[, colSums(tstats$tstatsfilt > 0) > 1 &
                      (.getOrigBaitNameFromComparison(colnames(tstats$tstatsfilt),
                                                      idmap = idmap) %in%
                           baitclass$Gene_name[baitclass$class500 != "Lost in 500 mM"])]
        tstat_onlybait <- tstat_onlybait[, colnames(tstat_onlybait) != "Untagged_plate"]
        tstat_several <- tstat_several[, colnames(tstat_several) != "Untagged_plate"]

        proteins_onlybait <-
            .getProteinNameFromSimplifiedComparison(colnames(tstat_onlybait))
        stopifnot(all(proteins_onlybait %in% rownames(tstat)))
        tstat_several <- tstat_several[!(rownames(tstat_several) %in% proteins_onlybait), ]

        ## Cluster proteins based on the correlation between their t-stat vectors
        dst <- sqrt(2 * (1 - cor(t(tstat_several), use = "pairwise.complete")))
        dst[is.na(dst)] <- sqrt(2 * (1 - (-1)))
        clst <- hclust(as.dist(dst))

        ## Reorder proteins by the clustering
        protorder <- c(rownames(tstat_several)[clst$order], proteins_onlybait)
        tstat <- tstat[protorder, ]

        ## Reorder conditions to keep a "diagonal" in the heatmap
        tstat <- reorderColsByRows(tstat)

        colsplit <- factor(
            structure(baitclass$class500[
                match(.getOrigBaitNameFromSimplifiedComparison(colnames(tstat),
                                                               idmap = idmap),
                      baitclass$Gene_name)],
                names = colnames(tstat)),
            levels = c("Retained in 500 mM", "Lost in 500 mM", ""))
        levels(colsplit) <- c("Interactions in\n500 mM NaCl IP-MS",
                              "No interactions in\n500 mM NaCl IP-MS", "")
    } else {
        stop("Unknown concentration")
    }
    return(list(tstat = tstat, colsplit = colsplit, nbrIntMat = nbrIntMat))
}

makeBaitHeatmapData <- function(sce, adjpthr, log2fcthr, nbrIntMat) {
    tc <- .getTestCols(sce, adjp_cutoff = adjpthr, logfc_cutoff = log2fcthr)$tstat
    colnames(tc) <- .getSimplifiedComparison(colnames(tc))
    shared_names <- intersect(rownames(tc),
                              .getProteinNameFromSimplifiedComparison(colnames(tc)))
    shared_names <- setdiff(shared_names, .getProteinFromOrigBait(nobait_150, idmap = idmap))
    stopifnot(all(shared_names %in% nbrIntMat$bait))
    nbrIntMat$complex <- NA
    for (cplx in c("MBF transcription complex", "CCAAT-binding factor complex", "Atf1-Pcr1")) {
        nbrIntMat$complex[.getPomBaseIdFromProtein(nbrIntMat$bait,
                                                   idmap = idmap) %in%
                              complexes[[cplx]]$pombase_ids] <- cplx
    }
    shared_names <- nbrIntMat[match(shared_names, nbrIntMat$bait), ] |>
        arrange(dbdfam, complex) |>
        pull(bait)
    mat <- tc[shared_names, match(shared_names, .getProteinNameFromSimplifiedComparison(colnames(tc)))]
    rownames(nbrIntMat) <- .getSimplifiedComparison(nbrIntMat$experiment)
    list(mat = mat, nbrIntMat = nbrIntMat)
}
