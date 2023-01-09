# SAM File Format
# ===============

module SAM

using BioGenerics

import BioAlignments
import BioGenerics: BioGenerics, isfilled, header
import BioGenerics.Exceptions: missingerror
import BioGenerics.Automa: State
import BioSequences
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..XAM: flag, XAMRecord, XAMReader, XAMWriter,
    ismapped, isprimary, ispositivestrand, isnextmapped #TODO: Deprecate import of flag queries. These were imported to preseve existing API.

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


include("metainfo.jl")
include("record.jl")
include("header.jl")
include("reader.jl")
include("readrecord.jl")
include("writer.jl")

end
