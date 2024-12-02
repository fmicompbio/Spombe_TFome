# obtain peak coordinates for a single transcription factor
# using publicly available data
#
# M. Stadler, 2024-11-08
#
# run this script in the same folder containing the repository from
# https://github.com/fmicompbio/Spombe_TFome


# load packages
library(GenomicRanges)
library(rtracklayer)

# load fused, filtered and annotated peak table
#    and convert to GRanges object
peaks <- as(read.csv("data/fused_peaks_filtered.csv.gz", row.names = 1), "GRanges")
table(peaks$peaktype)

# extract information about peak enrichments
isEnriched <- mcols(peaks)[, grep("^is_enr_in.", colnames(mcols(peaks)))]
colnames(isEnriched) <- sub("^is_enr_in.", "", colnames(isEnriched))
dim(isEnriched)

# print a list of all TF names (remove "Untagged")
tfNames <- setdiff(colnames(isEnriched), "Untagged")
tfNames

# iterate over TFs:
# for each TF, obtain peaks that are not "common peaks (ubiquitous)"
#     and output as bed file into the current working directory
for (tfSelect in tfNames) {
    peaksSelect <- peaks[isEnriched[, tfSelect] & peaks$peaktype != "common peaks (ubiquitous)"]
    message(tfSelect, " has ", length(peaksSelect), " specific/frequent peaks")
    if (length(peaksSelect) > 0) {
        tfBedfile <- paste0("peaks_", tfSelect, ".bed")
        export.bed(peaksSelect, tfBedfile)
    }
}
