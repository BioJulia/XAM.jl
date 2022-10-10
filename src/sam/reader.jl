# SAM Reader
# =========

mutable struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    header::Header
end

function Reader(state::State{S}) where {S <: TranscodingStream}

    rdr = Reader(state, Header())

    cs, ln = readheader!(rdr.state.stream, rdr.header, (sam_machine_header.start_state, rdr.state.linenum))

    rdr.state.state = sam_machine_body.start_state # Get the reader ready to read the body.
    rdr.state.linenum = ln
    rdr.state.filled = false

    if cs != -1 && cs != 0 #Note: the header is finished when the state machine fails to transition after a new line (state 1).
        throw(ArgumentError("Malformed SAM file header at line $(ln). Machine failed to transition from state $(cs)."))
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
    cs = index!(stream, record)
    if cs != 0
        throw(ArgumentError("Invalid SAM metadata. Machine failed to transition from state $(cs)."))
    end
    return record
end

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs = index!(stream, record)
    if cs != 0
        throw(ArgumentError("Invalid SAM record. Machine failed to transition from state $(cs)."))
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
function Base.read!(rdr::Reader, record::Record)

    cs, ln, found = readrecord!(rdr.state.stream, record, (rdr.state.state, rdr.state.linenum))

    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = found

    if found
        return record
    end

    if cs == 0 || eof(rdr.state.stream)
        throw(EOFError())
    end

    throw(ArgumentError("Malformed SAM file record at line $(ln). Machine failed to transition from state $(cs)."))

end
