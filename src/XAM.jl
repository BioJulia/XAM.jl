module XAM

export
    SAM,
    BAM

"""
    flag(record::Union{SAM.Record, BAM.Record})::UInt16

Get the bitwise flags of `record`. The returned value is a `UInt16` of each flag
being OR'd together. The possible flags are:

    0x0001 template having multiple segments in sequencing
    0x0002 each segment properly aligned according to the aligner
    0x0004 segment unmapped
    0x0008 next segment in the template unmapped
    0x0010 SEQ being reverse complemented
    0x0020 SEQ of the next segment in the template being reverse complemented
    0x0040 the first segment in the template
    0x0080 the last segment in the template
    0x0100 secondary alignment
    0x0200 not passing filters, such as platform/vendor quality controls
    0x0400 PCR or optical duplicate
    0x0800 supplementary alignment
"""
function flag end

include("sam/sam.jl")
include("bam/bam.jl")

using .SAM
using .BAM

end # module
