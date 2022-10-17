# BAM Reader
# ==========

"""
    BAM.Reader(input::IO; index=nothing)

Create a data reader of the BAM file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *bai* is supported) or BAI object
"""
mutable struct Reader{I<:Union{Nothing, BAI}} <: BioGenerics.IO.AbstractReader
    stream::BGZFStreams.BGZFStream
    header::SAM.Header
    start_offset::BGZFStreams.VirtualOffset
    refseqnames::Vector{String}
    refseqlens::Vector{Int}
    index::I
end

function Base.eltype(::Type{R}) where {R<:Reader}
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.stream
end

function Reader(input::IO; index=nothing)
    reader = init_bam_reader(input, index)
    return reader
end

function Base.show(io::IO, reader::Reader)
    println(io, summary(reader), ":")
      print(io, "  number of contigs: ", length(reader.refseqnames))
end

"""
    header(reader::Reader; fillSQ::Bool=false)::SAM.Header

Get the header of `reader`.

If `fillSQ` is `true`, this function fills missing "SQ" metainfo in the header.
"""
function header(reader::Reader; fillSQ::Bool=false)::SAM.Header
    header = reader.header
    if fillSQ
        if !isempty(findall(reader.header, "SQ"))
            throw(ArgumentError("SAM header already has SQ records"))
        end
        header = copy(header)
        for (name, len) in zip(reader.refseqnames, reader.refseqlens)
            push!(header, SAM.MetaInfo("SQ", ["SN" => name, "LN" => len]))
        end
    end
    return header
end

function Base.seek(reader::Reader, voffset::BGZFStreams.VirtualOffset)
    seek(reader.stream, voffset)
end

function Base.seekstart(reader::Reader)
    seek(reader.stream, reader.start_offset)
end

function Base.iterate(reader::Reader, nextone = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return copy(nextone), empty!(nextone)
end

# Initialize a BAM reader by reading the header section.
function init_bam_reader(input::BGZFStreams.BGZFStream, index = nothing)
    # magic bytes
    B = read(input, UInt8)
    A = read(input, UInt8)
    M = read(input, UInt8)
    x = read(input, UInt8)

    if B != UInt8('B') || A != UInt8('A') || M != UInt8('M') || x != 0x01
        error("input was not a valid BAM file")
    end

    # SAM header
    textlen = read(input, Int32)
    samreader = SAM.Reader(IOBuffer(read(input, textlen)))

    # reference sequences
    n_refs = read(input, Int32)
    refseqnames = Vector{String}(undef, n_refs)
    refseqlens = Vector{Int}(undef, n_refs)
    @inbounds for i in 1:n_refs
        namelen = read(input, Int32)
        data = read(input, namelen)
        seqname = unsafe_string(pointer(data))
        seqlen = read(input, Int32)
        refseqnames[i] = seqname
        refseqlens[i] = seqlen
    end

    voffset = isa(input.io, Base.AbstractPipe) ?
        BGZFStreams.VirtualOffset(0, 0) :
        BGZFStreams.virtualoffset(input)

    return Reader(
        input,
        samreader.header,
        voffset,
        refseqnames,
        refseqlens,
        init_bam_index(index)
    )
end

function init_bam_reader(input::IO, index = nothing)
    return init_bam_reader(BGZFStreams.BGZFStream(input), index)
end

init_bam_index(index::AbstractString) = BAI(index)
init_bam_index(index::BAI) = index
init_bam_index(index::Nothing) = nothing
init_bam_index(index) = error("unrecognizable index argument")

function _read!(reader::Reader, record)
    unsafe_read(
        reader.stream,
        pointer_from_objref(record),
        FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    if length(record.data) < dsize
        resize!(record.data, dsize)
    end
    unsafe_read(reader.stream, pointer(record.data), dsize)
    record.reader = reader
    return record
end
