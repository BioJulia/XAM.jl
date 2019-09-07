# SAM Reader
# =========

mutable struct Reader <: BioGenerics.IO.AbstractReader
    state::State
    header::Header

    function Reader(input::BufferedStreams.BufferedInputStream)
        reader = new(State(sam_header_machine.start_state, input), Header())
        readheader!(reader)
        reader.state.cs = sam_body_machine.start_state
        return reader
    end
end

"""
    SAM.Reader(input::IO)

Create a data reader of the SAM file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    return Reader(BufferedStreams.BufferedInputStream(input))
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

function Base.eltype(::Type{Reader})
    return Record
end
