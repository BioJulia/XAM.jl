# BAM File Format
# ===============

module BAM

using BioGenerics
using GenomicFeatures
using XAM.SAM
import ..XAM: flags, XAMRecord, XAMReader, XAMWriter,
	ismapped, isprimaryalignment, ispositivestrand, isnextmapped #TODO: Deprecate import of flag queries. These were imported to preseve existing API.

import BGZFLib
import MemoryViews
import BioAlignments
import Indexes
import BioSequences
import BioGenerics: isfilled, header
import BufferIO

import GenomicFeatures: eachoverlap

# Indexes.Chunk.start/stop are BGZFStreams.VirtualOffset (a 64-bit primitive type).
# Convert to BGZFLib.VirtualOffset without importing BGZFStreams.
@inline function _to_virtual_offset(vo)
    u = reinterpret(UInt64, vo)
    BGZFLib.VirtualOffset(u >> 16, u & 0xffff)
end

@inline function load_le(
            io::BufferIO.AbstractBufReader,
            ::Type{T}
    ) where T <: Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    # Fast path: Buffer has the data immediately available
    buffer = BufferIO.get_buffer(io)
    return if length(buffer) >= sizeof(T)
        value = GC.@preserve buffer htol(unsafe_load(Ptr{T}(pointer(buffer))))
        @inbounds BufferIO.consume(io, sizeof(T))
        value
    else
        load_le_slowpath(io, T)
    end
end

@noinline function load_le_slowpath(io::BufferIO.AbstractBufReader, T)
    # Second fastest path: A single fill_buffer (called via get_nonempty_buffer)
    # provides enough data, and then we load it immediately
    buffer = BufferIO.get_nonempty_buffer(io)
    sz = sizeof(T)
    buffer === nothing && throw(EOFError())
    if length(buffer) >= sz
        value = GC.@preserve buffer htol(unsafe_load(Ptr{T}(pointer(buffer))))
        @inbounds BufferIO.consume(io, sizeof(T))
        return value
    end

    # Very unlikely path: The T does not fit in the buffer of the IO, so
    # We need to store a temporary buffer.
    memory = MemoryViews.MemoryView(Memory{UInt8}(undef, sizeof(T)))
    remaining = memory
    while true
        buffer = @inbounds buffer[1:min(length(buffer), length(remaining))]
        @inbounds copyto!(remaining, buffer)
        @inbounds BufferIO.consume(io, length(buffer))
        remaining = @inbounds remaining[length(buffer) + 1 : end]
        if isempty(remaining)
            value = GC.@preserve memory htol(unsafe_load(Ptr{T}(pointer(memory))))
            return value
        end
        buffer = let
            b = BufferIO.get_nonempty_buffer(io)
            b === nothing && throw(EOFError())
            b
        end
    end
end

@inline function store_le(
        io::Union{BGZFLib.BGZFWriter, BGZFLib.SyncBGZFWriter},
        x::T
    ) where T <: Union{Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64}
    x = htol(x)
    # Fast path: Buffer already has room to store the T
    buffer = BufferIO.get_buffer(io)
    if length(buffer) >= sizeof(T)
        GC.@preserve buffer unsafe_store!(Ptr{T}(pointer(buffer)), x)
        @inbounds BufferIO.consume(io, sizeof(T))
        return nothing
    else
        store_le_slowpath(io, x)
    end
end

@noinline function store_le_slowpath(io::Union{BGZFLib.BGZFWriter, BGZFLib.SyncBGZFWriter}, x::T) where T
    BufferIO.grow_buffer(io)
    buffer = BufferIO.get_buffer(io)
    # BGZFLib documents that grow_buffer will do a shallow flush,
    # and that the buffer size is 2^16 bytes.
    # Some bytes are used for overhead, but this means we can guarantee around 2^16
    # bytes are available now.
    if length(buffer) < sizeof(T)
        error("Invalid BGZF implementation, too small buffer size after grow_buffer is called")
    end
    GC.@preserve buffer unsafe_store!(Ptr{T}(pointer(buffer)), x)
    @inbounds BufferIO.consume(io, sizeof(T))
    return nothing
end


include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
