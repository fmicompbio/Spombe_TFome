.capitalize <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

## Get the protein name (einprotID) corresponding to a given comparison
.getProteinNameFromComparison <- function(x) {
    ## x should be of the form sre1_150_plate_vs_highsaltcompl_sre1_150_plate.P.Value
    ## returns sre1
    .capitalize(sub("([^_]*)_[^_]*_([^_]*)_.*", "\\1", x))
}

## Get the bait name (original ID) corresponding to a given comparison
.getOrigBaitNameFromComparison <- function(x, idmap) {
    ## x should be of the form sre1_150_plate_vs_highsaltcompl_sre1_150_plate.P.Value
    ## returns the original bait/experiment name for sre1
    i <- sub("([^_]*)_[^_]*_([^_]*)_.*", "\\1", x)
    ifelse(.capitalize(i) == "Untagged", "Untagged",
           .capitalize(idmap$bait[match(i, .capitalize(idmap$unique_einprot_id))]))
}

## Simplify the comparison ID to just bait_{plate|tube}
.getSimplifiedComparison <- function(x) {
    ## x should be of the form sre1_150_plate_vs_highsaltcompl_sre1_150_plate.P.Value
    ## returns sre1_plate
    sub("([^_]*)_[^_]*_([^_]*)_.*", paste0("\\1", "_", "\\2"), x)
}

## Get the PomBase ID corresponding to a given comparison
.getPomBaseIdFromComparison <- function(x, idmap) {
    ## x should be of the form sre1_150_plate_vs_highsaltcompl_sre1_150_plate.P.Value
    ## returns SPBC19C2.09
    i <- sub("([^_]*)_[^_]*_([^_]*)_.*", "\\1", x)
    ifelse(.capitalize(i) == "Untagged", "Untagged",
           idmap$gene_stable_id[match(i, .capitalize(idmap$unique_einprot_id))])
}

## Get the UniProt ID corresponding to a given comparison
.getUniProtIdFromComparison <- function(x, idmap) {
    ## x should be of the form sre1_150_plate_vs_highsaltcompl_sre1_150_plate.P.Value
    ## returns SPBC19C2.09
    i <- sub("([^_]*)_[^_]*_([^_]*)_.*", "\\1", x)
    ifelse(.capitalize(i) == "Untagged", "Untagged",
           idmap$pombase_uniprot_id[match(i, .capitalize(idmap$unique_einprot_id))])
}

## Get the protein name (einprotID) corresponding to a given comparison
.getProteinNameFromSimplifiedComparison <- function(x) {
    ## x should be of the form sre1_plate
    ## returns sre1
    sub("([^_]*)_([^_]*)", "\\1", x)
}

## Get the bait name (original ID) corresponding to a given comparison
.getOrigBaitNameFromSimplifiedComparison <- function(x, idmap) {
    ## x should be of the form sre1_plate
    ## returns the original bait/experiment name for sre1
    i <- sub("([^_]*)_([^_]*)", "\\1", x)
    ifelse(.capitalize(i) == "Untagged", "Untagged",
           .capitalize(idmap$bait[match(i, .capitalize(idmap$unique_einprot_id))]))
}

## Get the PomBase ID corresponding to a given protein (unique einprot id)
.getPomBaseIdFromProtein <- function(x, idmap) {
    ## x should be of the form sre1
    ## returns SPBC19C2.09
    ifelse(.capitalize(x) == "Untagged", "Untagged",
           idmap$gene_stable_id[match(x, .capitalize(idmap$unique_einprot_id))])
}

## Get protein ID for a given original bait name
.getProteinFromOrigBait <- function(x, idmap) {
    ## x should be of the form sre1
    ## returns the corresponding protein name
    .capitalize(idmap$unique_einprot_id[match(x, .capitalize(idmap$bait))])
}
