module XAM

using BioGenerics
import BioGenerics: isfilled #Note: used by `ismapped`.

export
    SAM,
    BAM

abstract type XAMRecord end
abstract type XAMReader <: BioGenerics.IO.AbstractReader end
abstract type XAMWriter <: BioGenerics.IO.AbstractWriter end

include("flags.jl")

include("sam/sam.jl")
include("bam/bam.jl")
include("convert.jl")

using .SAM
using .BAM

end # module
