module XAM

export
    SAM,
    BAM

include("common.jl")
include("sam/sam.jl")
include("bam/bam.jl")

using .SAM
using .BAM

end # module
