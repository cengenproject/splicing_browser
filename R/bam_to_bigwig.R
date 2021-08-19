# Draft for visualization of raw alignments

library(GenomicAlignments)

awa <- coverage(readGAlignmentPairs("AWAr19.bam"))
nsm <- coverage(readGAlignmentPairs("NSMr59.bam"))

xx <- pmax(awa, nsm)
pmin(awa, nsm)
(awa + nsm)/2

pryr::object_size(nsm)
# 24,109,448 B  -> the Rle objects are small

rtracklayer::export.bw(awa, "export_awa.bw")
rtracklayer::export.bw(nsm, "export_NSM.bw")
rtracklayer::export.bw(pmax(awa,nsm), "export_max.bw")
rtracklayer::export.bw(pmin(awa,nsm), "export_min.bw")
rtracklayer::export.bw((awa+nsm)/2, "export_mean.bw")

