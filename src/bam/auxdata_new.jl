module NewAUXData

using ..BAM: BAM # Record, auxdata_position, data

# These are the numerical types supported by the BAM format.
const AUX_NUMBER_TYPES = Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32}

# We use a custom type to represent this for convenience, and also for
# performance. A string would need a heap allocation.
struct AUXTag
    x::Tuple{UInt8, UInt8}

    function AUXTag(x::UInt8, y::UInt8)
        is_valid_auxtag(x, y) || parse_tag_error()
        new((x, y))
    end
end

@noinline parse_tag_error() =
    throw(ArgumentError("Tags must conform to r\"[A-Za-z][A-Za-z0-9]\""))

# Per BAM specs, tags must conform to #[A-Za-z][A-Za-z0-9]
function is_valid_auxtag(x::UInt8, y::UInt8)
    digit = UInt8('0'):UInt8('9')
    upper = UInt8('A'):UInt8('Z')
    lower = UInt8('a'):UInt8('z')
    (in(x, upper) | in(x, lower)) & (in(y, digit) | in(y, upper) | in(y, lower))
end

AUXTag(x::Tuple{UInt8, UInt8}) = AUXTag(x...)

# Convert a char to the corresponding UInt8 if ASCII, else to 0xff
# This is more efficient than the Base implementation.
function to_uint8(x::Char)::UInt8
    u = reinterpret(UInt32, x)
    iszero(u & 0x80ffffff) ? 0xff : unsafe_u8(x)
end

function AUXTag(x::Char, y::Char)
    # If the chars are not bytes, to_uint8 will result in 0xff,
    # making the constructor error. 
    AUXTag(to_uint8(x), to_uint8(y))
end

function AUXTag(x::Union{String, SubString{String}})
    ncodeunits(x) == 2 || error("String must have two ASCII symbols")
    cu = codeunits(x)
    AUXTag(@inbounds cu[1], @inbounds cu[2])
end

# We don't want convert to allocate (since it's called implicitly), so
# we only accept string types that can efficiently be converted to tags
Base.convert(::Type{AUXTag}, x::Union{String, SubString{String}}) = AUXTag(x)

Base.show(io::IO, x::AUXTag) = write(io, "AUXTag(\"", x.x[1], x.x[2], "\")")

"""
    Hex(v::AbstractVector{UInt8})

Wrapper type around a byte vector. When this type is assigned to an AUXData,
it is encoded as a hex value (type tag `H`) as opposed to a raw byte array (
type tag `B`).
"""
struct Hex{V <: AbstractVector{UInt8}}
    x::V
end

struct AUXData <: AbstractDict{AUXTag, Any}
    # Same vector as the record's
    x::Vector{UInt8}
    # First index of the AUX data in the record
    start::Int
end

AUXData() = AUXData(UInt8[], 1)

function AUXData(itr)
    y = AUXData()
    for (k, v) in itr
        tag = convert(AUXTag, k)
        setindex!(y, v, tag)
    end
    y
end

AUXData(x::BAM.Record) = AUXData(x.data, BAM.auxdata_position(x))

# For this constructor, we can rely on the keys being unique, so we can
# use the more efficient setindex_nonexisting!
function AUXData(d::AbstractDict)
    y = AUXData()
    for (k, v) in pairs(d)
        tag = convert(AUXTag, k)
        setindex_nonexisting!(y, v, tag)
    end
    y
end

# Because checking the size is O(N).
Base.IteratorSize(::Type{AUXData}) = Base.SizeUnknown()
Base.length(x::AUXData) = sum(1 for i in x; init=0)

@noinline invalid_aux() = error("Invalid AuxData")

# Step function to get the index of the next tag in the vector.
function get_next_tag(x::AUXData, i::Int)::Union{Nothing, Tuple{AUXTag, Int}}
    v = x.x
    i > length(v) && return nothing
    # Two bytes tag plus 1 byte type tag
    data_length = length(v) - (i + 2)
    # If there is zero bytes for any values. We error here to avoid reading OOB
    # for the tag and type tag bytes
    data_length < 1 && invalid_aux() # TODO: Error values instead of throwing?
    tag = AUXTag(@inbounds v[i], @inbounds v[i + 1])
    type_tag = @inbounds v[i + 2]
    start = i + 3

    # The following could be done with dispatch. However, this would lead to
    # type instability since the dispatch would depend on the value of the
    # type tag byte. So, if/else spam is (much) more efficient.

    # One byte values
    if type_tag in (UInt8('C'), UInt8('c'), UInt8('A'))
        return (tag, start + 1)
        # Two byte values
    elseif type_tag in (UInt8('S'), UInt8('s'))
        data_length < 2 && invalid_aux()
        return (tag, start + 2)
        # Four byte values
    elseif type_tag in (UInt8('I'), UInt8('i'), UInt8('f'))
        data_length < 4 && invalid_aux()
        return (tag, start + 4)
        # Null-terminated values
    elseif type_tag in (UInt8('Z'), UInt8('H'))
        zeropos = findnext(iszero, v, start)
        isnothing(zeropos) && invalid_aux()
        return (tag, zeropos + 1)
        # Arrays
    elseif type_tag == UInt8('B')
        # Minimum data length for empty array:
        # Array element type byte plus 4 for array length
        data_length < 5 && invalid_aux()
        eltype_tag = v[i + 3]
        # Don't use dispatch here for efficiency to avoid type instability
        eltype_size = if eltype_tag in (UInt8('C'), UInt8('c'))
            1
        elseif eltype_tag in (UInt8('S'), UInt8('s'))
            2
        elseif eltype_tag in (UInt8('I'), UInt8('i'), UInt8('f'))
            4
        else
            invalid_aux()
        end
        # Note: All BAM integers are little endian so we can do this
        n_elements = @inbounds begin
            v[i + 4] % UInt32 |
            (v[i + 5] % UInt32) << 8 |
            (v[i + 6] % UInt32) << 16 |
            (v[i + 7] % UInt32) << 24
        end
        len = n_elements * eltype_size
        data_length < len + 5 && invalid_aux()
        return (tag, start + 5 + len)
    else
        invalid_aux()
    end
end

const ELTYPE_DICT = Dict(
    UInt8('C') => UInt8,
    UInt8('c') => Int8,
    UInt8('S') => UInt16,
    UInt8('s') => Int16,
    UInt8('I') => UInt32,
    UInt8('i') => Int32,
    UInt8('f') => Float32,
)

# BAM files only allow these bytes for chars in AUXData
is_printable_char(x::UInt8) = in(x, UInt8('!'):UInt8('~'))
function is_printable_char(c::Char)
    u = reinterpret(UInt32, c)
    iszero(u & 0x80ffffff) & is_printable_char(unsafe_u8(c))
end
unsafe_u8(c::Char) = (reinterpret(UInt32, c) >> 24) % UInt8

# Get the value from AUXData using the span.
# The span MUST be a valid span, so this function is to be used for internal use only.
function materialize(x::AUXData, span::UnitRange{Int})
    v = x.x
    # Skip first two bytes: The AUXTag
    i = first(span) + 2
    type_tag = v[i]
    i += 1
    if type_tag == UInt8('C')
        @inbounds v[i]
    elseif type_tag == UInt8('c')
        v[i] % Int8
    elseif type_tag == UInt8('A')
        c = v[i]
        is_printable_char(c) || invalid_aux()
        reinterpret(Char, (c % UInt32) << 24)
    else
        GC.@preserve v begin
            if type_tag == UInt8('s')
                ltoh(unsafe_load(Ptr{Int16}(pointer(v, i))))
            elseif type_tag == UInt8('S')
                ltoh(unsafe_load(Ptr{UInt16}(pointer(v, i))))
            elseif type_tag == UInt8('i')
                ltoh(unsafe_load(Ptr{Int32}(pointer(v, i))))
            elseif type_tag == UInt8('I')
                ltoh(unsafe_load(Ptr{UInt32}(pointer(v, i))))
            elseif type_tag == UInt8('f')
                ltoh(unsafe_load(Ptr{Float32}(pointer(v, i))))
            elseif type_tag == UInt8('Z')
                # Compensate for null terminator byte
                validate_aux_value(String(v[i:(last(span) - 1)]))
            elseif type_tag == UInt8('H')
                # Compensate for null terminator byte
                hex2bytes(@view v[i:(last(span) - 1)])
            elseif type_tag == UInt8('B')
                T = ELTYPE_DICT[v[i]]
                collect(reinterpret(T, v[(i + 5):last(span)]))
            else
                error() # unreachable!
            end
        end
    end
end

# Get the span of a key if it is present, or nothing if not. Has O(N) time.
function get_span(x::AUXData, k)
    tag = convert(AUXTag, k)
    start = x.start
    while true
        nexttag = get_next_tag(x, start)
        nexttag === nothing && return nothing
        (seen_tag, next_index) = nexttag
        seen_tag == tag && return start:(next_index - 1)
        start = next_index
    end
end

function Base.delete!(x::AUXData, key)
    tag = convert(AUXTag, key)
    ind = get_span(x, tag)
    if !isnothing(ind)
        deleteat!(x.x, ind)
    end
    x
end

Base.haskey(x::AUXData, k) = !isnothing(get_span(x, k))

function Base.get(x::AUXData, k, default)
    ind = get_span(x, k)
    isnothing(ind) ? default : materialize(x, ind)
end

function Base.iterate(x::AUXData, state::Int=Int(x.start))
    # This first call is always safe and will return nothing if x is empty.
    nexttag = get_next_tag(x, state)
    isnothing(nexttag) && return nothing
    (tag, next_index) = nexttag
    value = materialize(x, state:(next_index - 1))
    (tag => value, next_index)
end

function Base.setindex!(x::AUXData, val, k)
    tag = convert(AUXTag, k)
    # Check if the tag already exists, and if so, where it is
    span = get_span(x, tag)
    # If it does not exist, insert it at the end
    if span === nothing
        setindex_nonexisting!(x, val, k)
    else
        # Else, if tag exists, check how many bytes is needed for the data.
        # We need to write the type tag and the 
        v = x.x
        value = as_aux_value(val)
        # Needed: bytes for data plus one for the type tag
        needed = 1 + bytes_needed(value)
        # Existing data: The entire span minus the two for the tag
        L = length(span) - 2
        # If the old value took up less space, we need to allocate extra space,
        # and then move all bytes after the current value to make room for the new one
        if L < needed
            N = length(v) - last(span)
            resize!(v, length(v) + needed - L)
            copyto!(v, last(span) + 1 + needed - L, v, last(span) + 1, N)
            # If the old value took up more space, we need to delete some data
            # and shrink the AUXData's stop value to match
        elseif L > needed
            n_remove = L - needed
            deleteat!(v, (first(span) + 2):(first(span) + 2 + n_remove - 1))
        end
        # Write in the new type tag, then the value itself
        v[first(span) + 2] = get_type_tag(typeof(value))
        write_data!(v, value, first(span) + 3)
    end
end

# Does not check for existing key
function setindex_nonexisting!(x::AUXData, v, k)
    tag = convert(AUXTag, k)
    # Convert the value to one of a type that can be stored in the AUXData,
    # and validate its content.
    val = validate_aux_value(as_aux_value(v))
    vec = x.x
    write_index = length(vec) + 1
    # Bytes needed is just for the data payload itself, so we need 2 for the tag
    # and one for the type tag
    resize!(vec, length(vec) + 2 + 1 + bytes_needed(val))
    # Write tag, then type tag, then data
    unsafe_write_tag!(vec, tag, write_index)
    vec[write_index + 2] = get_type_tag(typeof(val))
    write_data!(vec, val, write_index + 3)
    x
end

as_aux_type(T::Type{<:AUX_NUMBER_TYPES}) = T
as_aux_type(::Type{<:Real}) = Float32
as_aux_type(::Type{<:Integer}) = Int32

as_aux_value(x::Union{Char, String, SubString{String}, Hex, Vector{<:AUX_NUMBER_TYPES}}) = x
as_aux_value(s::Union{String, SubString{String}}) = s
as_aux_value(s::AbstractString) = String(s)::String

function as_aux_value(x::Real)
    T = as_aux_type(typeof(x))
    T(x)::T
end

function as_aux_value(v::AbstractVector{T}) where {T <: Real}
    E = as_aux_type(T)::AUX_NUMBER_TYPES
    Vector{E}(v)::Vector{E}
end

# Default implementation: Check nothing
validate_aux_value(x) = x

function validate_aux_value(c::Char)
    is_printable_char(c) || error() # TODO
    c
end

# Strings in BAM can be printable chars, or space, but nothing else.
function validate_aux_value(s::Union{String, SubString{String}})
    isvalid = true
    for i in codeunits(s)
        isvalid &= (i == UInt8(' ')) | is_printable_char(i)
    end
    isvalid || error() # TODO
    s
end

bytes_needed(x::Union{Int8, UInt8, Char}) = 1
bytes_needed(x::Union{Int16, UInt16}) = 2
bytes_needed(x::Union{Int32, UInt32, Float32}) = 4
bytes_needed(x::AbstractString) = ncodeunits(x) + 1 # null byte at end
bytes_needed(x::Hex) = 2 * length(x.x) + 1 # null byte

function bytes_needed(x::Vector{<:AUX_NUMBER_TYPES})
    # Element type tag + 4 bytes for length, plus data itself
    1 + 4 + sizeof(x)
end

function unsafe_write_tag!(v::Vector{UInt8}, tag::AUXTag, i::Integer)
    @inbounds v[i] = tag.x[1]
    @inbounds v[i + 1] = tag.x[2]
    v
end

# These are the type tags BAM use to determine the value of the serialized data
get_type_tag(::Type{UInt8}) = UInt8('C')
get_type_tag(::Type{Int8}) = UInt8('c')
get_type_tag(::Type{UInt16}) = UInt8('S')
get_type_tag(::Type{Int16}) = UInt8('s')
get_type_tag(::Type{UInt32}) = UInt8('I')
get_type_tag(::Type{Int32}) = UInt8('i')
get_type_tag(::Type{Float32}) = UInt8('f')
get_type_tag(::Type{Char}) = UInt8('A')
get_type_tag(::Type{<:Union{String, SubString{String}}}) = UInt8('Z')
get_type_tag(::Type{<:AbstractVector}) = UInt8('B')
get_type_tag(::Type{<:Hex}) = UInt8('H')

function write_data!(v::Vector{UInt8}, val::Char, i::Int)
    v[i] = unsafe_u8(val)
end

function write_data!(v::Vector{UInt8}, val::Union{Int8, UInt8}, i::Int)
    v[i] = reinterpret(UInt8, val)
end

function write_data!(v::Vector{UInt8}, val::AUX_NUMBER_TYPES, i::Int)
    GC.@preserve v begin
        unsafe_store!(Ptr{typeof(val)}(pointer(v, i)), htol(val))
    end
end

function write_data!(v::Vector{UInt8}, val::Union{String, SubString{String}}, i::Int)
    for (n, u) in enumerate(codeunits(val))
        v[i + n - 1] = u
    end
    # Null byte termination
    v[i + ncodeunits(val)] = 0x00
end

hexencode_nibble(u::UInt8) = u < 0x0a ? UInt8('0') + u : UInt8('A') - 0x0a + u
function write_data!(v::Vector{UInt8}, val::Hex, i::Int)
    for (byte_no, byte) in enumerate(val.x)
        v[2 * (byte_no - 1) + i] = hexencode_nibble(byte >> 4)
        v[2 * (byte_no - 1) + i + 1] = hexencode_nibble(byte & 0x0f)
    end
    # Null byte termination
    v[i + 2 * length(val.x)] = 0x00
end

function write_data!(v::Vector{UInt8}, val::Vector{E}, i::Int) where {E <: AUX_NUMBER_TYPES}
    # Write eltype tag first
    v[i] = get_type_tag(E)
    GC.@preserve v begin
        ptr = Ptr{E}(pointer(v, i + 1))
        # Write the 32-bit array length
        unsafe_store!(Ptr{UInt32}(ptr), htol(length(val) % UInt32))
        ptr += 4
        for e in val
            unsafe_store!(ptr, htol(E(e)::E))
            ptr += sizeof(E)
        end
    end
end

end # module
