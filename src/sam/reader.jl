# SAM Reader
# =========

mutable struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    header::Header
end

function Reader(state::State{S}) where {S <: TranscodingStream}

    rdr = Reader(state, Header())

    cs, ln, f = readheader!(rdr.state.stream, rdr.header, (sam_machine_header.start_state, rdr.state.linenum))

    rdr.state.state = sam_machine_body.start_state
    rdr.state.linenum = ln
    rdr.state.filled = false

    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("Malformed SAM file at line $(ln)."))
    end
    
    return rdr
end

"""
    SAM.Reader(input::IO)

Create a data reader of the SAM file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)

    if input isa TranscodingStream
        return Reader(State(input, 1, 1, false))
    end

    stream = TranscodingStreams.NoopStream(input)

    return Reader(State(stream, 1, 1, false))

end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

"""
    header(reader::Reader)::Header

Get the header of `reader`.
"""
function header(reader::Reader)::Header
    return reader.header
end

function Base.close(reader::Reader)
    if reader.state.stream isa IO
        close(reader.state.stream)
    end
    return nothing
end

function index!(record::MetaInfo)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    found = index!(stream, record)
    if !found
        throw(ArgumentError("invalid SAM metadata"))
    end
    return record
end

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    found = index!(stream, record)
    if !found
        throw(ArgumentError("invalid SAM record"))
    end
    return record
end

function Base.iterate(reader::Reader, nextone::Record = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return copy(nextone), empty!(nextone)
end

"""
    read!(rdr::Reader, rec::Record)

Read a `Record` into `rec`; overwriting or adding to existing field values.
It is assumed that `rec` is already initialized or empty.
"""
function Base.read!(rdr::Reader, rec::Record)

    cs, ln, f = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum))

    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f

    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("Malformed SAM file at line $(ln)."))
    end

    return rec
end
