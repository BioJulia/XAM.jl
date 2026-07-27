# BAM Writer
# ==========

"""
    BAM.Writer(output::T, header::SAM.Header)

Create a data writer of the BAM file format.

# Arguments
* `output`: data sink
* `header`: SAM header object
"""
mutable struct Writer{T} <: XAMWriter
    stream::T
end

function Writer(stream::T, header::SAM.Header) where T
    refseqnames = String[]
    refseqlens = Int[]
    for metainfo in findall(header, "SQ")
        push!(refseqnames, metainfo["SN"])
        push!(refseqlens, parse(Int, metainfo["LN"]))
    end
    write_header(stream, header, refseqnames, refseqlens)
    return Writer(stream)
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, record::Record)
    n = 0
    n += unsafe_write(writer.stream, pointer_from_objref(record), UInt64(FIXED_FIELDS_BYTES))
    n += unsafe_write(writer.stream, pointer(record.data), UInt64(data_size(record)))
    return n
end

function write_header(stream, header, refseqnames, refseqlens)
    @assert length(refseqnames) == length(refseqlens) "Lengths of refseq names and lengths must match."
    
    # magic bytes
    write(stream, "BAM\1")

    # SAM header
    buf = IOBuffer()
    l = write(SAM.Writer(buf), header)
    store_le(stream, Int32(l))
    write(stream, take!(buf))

    # reference sequences
    store_le(stream, Int32(length(refseqnames)))
    for (seqname, seqlen) in zip(refseqnames, refseqlens)
        namelen = length(seqname)
        store_le(stream, Int32(namelen + 1))
        write(stream, seqname, '\0')
        store_le(stream, Int32(seqlen))
    end

    return nothing
end
