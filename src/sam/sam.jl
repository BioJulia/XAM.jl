# SAM File Format
# ===============

module SAM

using BioGenerics

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioAlignments
import BioGenerics: BioGenerics, isfilled, header
import BioGenerics.Exceptions: missingerror
import BioGenerics.Automa: State
import BioSequences
import TranscodingStreams: TranscodingStreams, TranscodingStream

using Printf: @sprintf


#TODO: update import BioCore.RecordHelper: unsafe_parse_decimal
# r"[0-9]+" must match `data[range]`.
function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where {T<:Unsigned}
    x = zero(T)
    @inbounds for i in range
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return x
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where {T<:Signed}
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

# #TODO: update BioCore.Ragel.State (will likely change with TrnscodingStreams).
# import BufferedStreams: BufferedStreams, BufferedInputStream
# # A type keeping track of a ragel-based parser's state.
# mutable struct State{T<:BufferedInputStream}
#     stream::T      # input stream
#     cs::Int        # current DFA state of Ragel
#     linenum::Int   # line number: parser is responsible for updating this
#     finished::Bool # true if finished (regardless of where in the stream we are)
# end

# function State(initstate::Int, input::BufferedInputStream)
#     return State(input, initstate, 1, false)
# end


include("flags.jl")
include("metainfo.jl")
include("record.jl")
include("header.jl")
include("reader.jl")
include("readrecord.jl")
include("writer.jl")

end
