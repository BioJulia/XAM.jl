module XAM

export
    SAM,
    BAM,
    header,
    eachoverlap,
    isfilled,
    seqname,
    hasseqname,
    leftposition,
    hasleftposition,
    rightposition,
    hasrightposition,
    sequence,
    hassequence

import BioCore: BioCore, distance, header, isfilled, seqname, hasseqname, sequence, hassequence, leftposition, rightposition, hasleftposition, hasrightposition
import GenomicFeatures: eachoverlap

include("sam/sam.jl")
include("bam/bam.jl")

end # module
