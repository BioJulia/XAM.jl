module XAM

export
    SAM,
    BAM

include("sam/sam.jl")
include("bam/bam.jl")

using .SAM
using .BAM

end # module
