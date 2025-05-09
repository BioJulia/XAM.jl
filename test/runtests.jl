using Test
using Documenter

using BioGenerics
using FormatSpecimens
using GenomicFeatures
using XAM

import BioAlignments: Alignment, AlignmentAnchor, OP_START, OP_MATCH, OP_DELETE
import BGZFStreams: BGZFStream
import BioGenerics.Exceptions: MissingFieldException
import BioSequences: @dna_str, @aa_str


# Generate a random range within `range`.
function randrange(range)
    x = rand(range)
    y = rand(range)
    if x < y
        return x:y
    else
        return y:x
    end
end


@testset "XAM" begin
    include("test_sam.jl")
    include("test_bam.jl")
    include("test_issues.jl")
    include("test_crosscheck.jl")
    include("test_convert.jl")
    include("test_flags.jl")

    # Include doctests.
    DocMeta.setdocmeta!(XAM, :DocTestSetup, :(using XAM); recursive=true)
    doctest(XAM; manual = false)
end
