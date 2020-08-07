# BAM Record
# ==========

"""
    BAM.Record()

Create an unfilled BAM record.
"""
mutable struct Record
    # fixed-length fields (see BMA specs for the details)
    block_size::Int32
    refid::Int32
    pos::Int32
    l_read_name::UInt8
    mapq::UInt8
    bin::UInt16
    n_cigar_op::UInt16
    flag::UInt16
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32
    # variable length data
    data::Vector{UInt8}
    reader::Union{Reader, Nothing}

    function Record()
        # An empty record is a legitimate BAM record.
        block_size = FIXED_FIELDS_BYTES - sizeof(Int32) + 1
        # Only the null terminator of the query name sequence
        data = UInt8[0x00]
        return new(block_size, -1, -1, 0x01, 0xff, 0, 0, 0x0004, 0, -1, -1, 0, data, nothing)
    end
end

# the data size of fixed-length fields (block_size-tlen)
const FIXED_FIELDS_BYTES = 36

function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    length(data) < FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    record = Record()
    dst_pointer = Ptr{UInt8}(pointer_from_objref(record))
    unsafe_copyto!(dst_pointer, pointer(data), FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    resize!(record.data, dsize)
    length(data) < dsize + FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    unsafe_copyto!(record.data, 1, data, FIXED_FIELDS_BYTES + 1, dsize)
    return record
end

function Base.:(==)(a::Record, b::Record)
    return a.block_size == b.block_size &&
        a.refid         == b.refid &&
        a.pos           == b.pos &&
        a.l_read_name   == b.l_read_name &&
        a.mapq          == b.mapq &&
        a.bin           == b.bin &&
        a.n_cigar_op    == b.n_cigar_op &&
        a.flag          == b.flag &&
        a.l_seq         == b.l_seq &&
        a.next_refid    == b.next_refid &&
        a.next_pos      == b.next_pos &&
        a.tlen          == b.tlen &&
        a.data[1:data_size(a)] == b.data[1:data_size(b)]
end

function Base.copy(record::Record)
    copy = Record()
    GC.@preserve copy record begin
    	dst_pointer = Ptr{UInt8}(pointer_from_objref(copy))
    	src_pointer = Ptr{UInt8}(pointer_from_objref(record))
    	unsafe_copyto!(dst_pointer, src_pointer, FIXED_FIELDS_BYTES)
    end
    copy.data       = record.data[1:data_size(record)]
    copy.reader     = record.reader
    return copy
end

function Base.empty!(record::Record)
    block_size         =  FIXED_FIELDS_BYTES - sizeof(Int32) + 1
    record.block_size  = block_size
    record.refid       = -1
    record.pos         = -1
    record.l_read_name = 0x01
    record.mapq        = 0xff
    record.bin         = 0
    record.n_cigar_op  = 0
    record.flag        = 0x0004
    record.l_seq       = 0
    record.next_refid  = -1
    record.next_pos    = -1
    record.tlen        = 0
    record.data[1]     = 0x00

    #Note: data will be overwritten and indexed using data_size.
    return record
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    println(io)
    println(io, "      template name: ", tempname(record))
    println(io, "               flag: ", flag(record))
    println(io, "       reference ID: ", refid(record))
    println(io, "           position: ", position(record))
    println(io, "    mapping quality: ", mappingquality(record))
    println(io, "              CIGAR: ", cigar(record))
    println(io, "  next reference ID: ", nextrefid(record))
    println(io, "      next position: ", nextposition(record))
    println(io, "    template length: ", templength(record))
    println(io, "           sequence: ", sequence(record))
    # TODO: pretty print base quality
    println(io, "       base quality: ", quality(record))
      print(io, "     auxiliary data:")
    for field in keys(auxdata(record))
        print(io, ' ', field, '=', record[field])
    end
end

function Base.read!(reader::Reader, record::Record)
    return _read!(reader, record)
end


# Accessor Fuctions
# -----------------

"""
    flag(record::Record)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::Record)::UInt16
    return record.flag
end

"""
    ismapped(record::Record)::Bool

Test if `record` is mapped.
"""
function ismapped(record::Record)::Bool
    return flag(record) & SAM.FLAG_UNMAP == 0
end

"""
    isprimary(record::Record)::Bool

Test if `record` is a primary line of the read.

This is equivalent to `flag(record) & 0x900 == 0`.
"""
function isprimary(record::Record)::Bool
    return flag(record) & 0x900 == 0
end

"""
    ispositivestrand(record::Record)::Bool

Test if `record` is aligned to the positive strand.

This is equivalent to `flag(record) & 0x10 == 0`.
"""
function ispositivestrand(record::Record)::Bool
    flag(record) & 0x10 == 0
end

"""
    refid(record::Record)::Int

Get the reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is 0 for a record without a mapping position.

See also: `BAM.rname`
"""
function refid(record::Record)
    hasrefid(record) || return nothing
	return Int(record.refid + 1)
end

function hasrefid(record::Record)
	record.refid > -1
end

"""
    refname(record::Record)::String

Get the reference sequence name of `record`.

See also: `BAM.refid`
"""
function refname(record::Record)
    hasrefname(record) || return nothing
    return record.reader.refseqnames[refid(record)]
end

"""
    reflen(record::Record)::Int

Get the length of the reference sequence this record applies to.
"""
function reflen(record::Record)
	id = refid(record)
	if (id === nothing) | (record.reader === nothing)
		return nothing
	end
    return record.reader.refseqlens[id]
end

function hasrefname(record::Record)
    return hasrefid(record) & (record.reader !== nothing)
end

"""
    position(record::Record)::Int

Get the 1-based leftmost mapping position of `record`.
"""
function position(record::Record)
	return hasposition(record) ? record.pos + 1 : nothing
end

function hasposition(record::Record)
	return record.pos > -1
end

"""
    rightposition(record::Record)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::Record)
    pos = position(record)
    return pos === nothing ? nothing : pos + alignlength(record) - 1
end

function hasrightposition(record::Record)
    return record.pos > -1
end

"""
    isnextmapped(record::Record)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::Record)
    return flag(record) & SAM.FLAG_MUNMAP == 0
end

"""
    nextrefid(record::Record)::Int

Get the next/mate reference sequence ID of `record`.
"""
function nextrefid(record::Record)
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    if !ispaired || record.next_refid == -1
        return nothing
    end
    return Int(record.next_refid + 1)
end

function hasnextrefid(record::Record)
    return record.next_refid > -1
end

"""
    nextrefname(record::Record)::String

Get the reference name of the mate/next read of `record`.
"""
function nextrefname(record::Record)
    hasnextrefname(record) || return nothing
    id = record.next_refid
    return record.reader.refseqnames[id]
end

function hasnextrefname(record::Record)
    return (record.next_refid > -1) & (record.reader !== nothing)
end

"""
    nextposition(record::Record)::Int

Get the 1-based leftmost mapping position of the next/mate read of `record`.
"""
function nextposition(record::Record)
    hasnextposition(record) || return nothing
    return record.next_pos + 1
end

function hasnextposition(record::Record)
    return record.next_pos > -1
end

"""
    mappingquality(record::Record)::UInt8

Get the mapping quality of `record`.
"""
function mappingquality(record::Record)
    hasmappingquality(record) || return nothing
    return record.mapq
end

function hasmappingquality(record::Record)
    return record.mapq < 0xff
end

"""
    n_cigar_op(record::Record, checkCG::Bool = true)

Return the number of operations in the CIGAR string of `record`.

Note that in the BAM specification, the field called `cigar` typically stores the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long, and due to  some constraints of the BAM format, the actual cigar string is stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar string.

Calling this method with `checkCG` set to `true` (default) this method will always yield the number of operations in the true cigar string, because this is probably what you want, the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still want to get the number of operations in the `cigar` field of the BAM record, then set `checkCG` to `false`.
"""
function n_cigar_op(record::Record, checkCG::Bool = true)
    return cigar_position(record, checkCG)[2]
end

"""
    cigar(record::Record)::String

Get the CIGAR string of `record`.

Note that in the BAM specification, the field called `cigar` typically stores the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long, and due to  some constraints of the BAM format, the actual cigar string is stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar string.

Calling this method with `checkCG` set to `true` (default) this method will always yield the true cigar string, because this is probably what you want the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still want to access the pseudo-cigar that is stored in the `cigar` field of the BAM record, then you can set checkCG to `false`.

See also `BAM.cigar_rle`.
"""
function cigar(record::Record, checkCG::Bool = true)::String
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(record, checkCG)...)
        print(buf, len, convert(Char, op))
    end
    return String(take!(buf))
end

"""
    cigar_rle(record::Record, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}

Get a run-length encoded tuple `(ops, lens)` of the CIGAR string in `record`.

Note that in the BAM specification, the field called `cigar` typically stores the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long, and due to  some constraints of the BAM format, the actual cigar string is stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar string.

Calling this method with `checkCG` set to `true` (default) this method will always yield the true cigar string, because this is probably what you want the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still want to access the pseudo-cigar that is stored in the `cigar` field of the BAM record, then you can set checkCG to `false`.

See also `BAM.cigar`.
"""
function cigar_rle(record::Record, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}
    idx, nops = cigar_position(record, checkCG)
    ops, lens = extract_cigar_rle(record.data, idx, nops)
    return ops, lens
end

function extract_cigar_rle(data::Vector{UInt8}, index, n)
    ops = Vector{BioAlignments.Operation}(undef, n)
    lens = Vector{Int}(undef, n)
    GC.@preserve data @inbounds for i in 1:n
        x = unsafe_load(Ptr{UInt32}(pointer(data, index)))
        ops[i] = BioAlignments.Operation(x & 0x0F)
        lens[i] = x >>> 4
        index += 4
    end
    return ops, lens
end

function cigar_position(record::Record, checkCG::Bool = true)::Tuple{Int, Int}
    cigaridx, nops = seqname_length(record) + 2, record.n_cigar_op
    # If a read has more than 65k ops, it's encoded in an AUX field, and the ops
    # is set to kSmN. So first return if it's not that
    if !checkCG
        return cigaridx, nops
    end
    if nops != 2
        return cigaridx, nops
    end
    x = unsafe_load(Ptr{UInt32}(pointer(record.data, cigaridx)))
    if x != UInt32(seqlength(record) << 4 | 4)
        return cigaridx, nops
    end
    # Else we go fetch it from the AUX fields.
    start = auxdata_position(record)
    stop = data_size(record)
    tagidx = findauxtag(record.data, start, stop, UInt8('C'), UInt8('G'))
    if tagidx == 0
        return cigaridx, nops
    end
    # Tag exists, validate type is BI.
    typ = unsafe_load(Ptr{UInt16}(pointer(record.data, tagidx += 2)))
    if typ != (UInt16('I') << 8 | UInt16('B'))
        return cigaridx, nops
    end
    # If got this far, the CG tag is valid and contains the cigar.
    # Get the true n_cigar_ops, and return it and the idx of the first
    nops = UInt32(unsafe_load(Ptr{Int32}(pointer(record.data, tagidx += 2))))
    tagidx += 4
    return tagidx, nops
end

"""
    alignment(record::Record)::BioAlignments.Alignment

Get the alignment of `record`.
"""
function alignment(record::Record)::BioAlignments.Alignment
    if !ismapped(record)
        return BioAlignments.Alignment(BioAlignments.AlignmentAnchor[])
    end
    seqpos = 0
    refpos = position(record) - 1
    anchors = [BioAlignments.AlignmentAnchor(seqpos, refpos, BioAlignments.OP_START)]
    for (op, len) in zip(cigar_rle(record)...)
        if BioAlignments.ismatchop(op)
            seqpos += len
            refpos += len
        elseif BioAlignments.isinsertop(op)
            seqpos += len
        elseif BioAlignments.isdeleteop(op)
            refpos += len
        else
            error("operation $(op) is not supported")
        end
        push!(anchors, BioAlignments.AlignmentAnchor(seqpos, refpos, op))
    end
    return BioAlignments.Alignment(anchors)
end

function hasalignment(record::Record)
    return ismapped(record)
end

"""
    alignlength(record::Record)::Int

Get the alignment length of `record`.
"""
function alignlength(record::Record)
    idx, nops = cigar_position(record, false)
    length = 0
    data = record.data
    GC.@preserve data for i in 1:nops
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, idx)), i)
        op = BioAlignments.Operation(x & 0x0F)
        if BioAlignments.ismatchop(op) | BioAlignments.isdeleteop(op)
            length += x >> 4
        end
    end
    return length
end

"""
    tempname(record::Record)::String

Get the query template name of `record`.
"""
function tempname(record::Record)
    hastempname(record) || return nothing
    return unsafe_string(pointer(record.data), seqname_length(record))
end

function hastempname(record::Record)
    return seqname_length(record) > 0
end

"""
    templength(record::Record)::Int

Get the template length of `record`.
"""
function templength(record::Record)
    hastemplength(record) || return nothing
    return record.tlen
end

function hastemplength(record::Record)
    return record.tlen > 0
end

"""
    sequence(record::Record)::BioSequences.LongDNASeq

Get the segment sequence of `record`.
"""
function sequence(record::Record)
    hassequence(record) || return nothing
    seqlen = seqlength(record)
    data = Vector{UInt64}(undef, cld(seqlen, 16))
    index = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + 1
    src::Ptr{UInt64} = pointer(record.data, index)
    for i in 1:lastindex(data)
        # copy data flipping high and low nybble
        x = unsafe_load(src, i)
        data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
    end
    return BioSequences.LongDNASeq(data, 1:seqlen, false)
end

function hassequence(record::Record)
    return hasseqlength(record)
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)
    hasseqlength(record) || return nothing
    return record.l_seq % Int
end

function hasseqlength(record::Record)
    return !iszero(record.l_seq)
end

"""
    quality(record::Record)

Get the base quality of `record`.
"""
function quality(record::Record)
    hasseqlength(record) || return nothing
    seqlen = seqlength(record)
    offset = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + cld(seqlen, 2)
    quals = record.data[(1+offset):(seqlen+offset)]
    all(i == 0xff for i in quals) && return nothing
    return quals
end

function hasquality(record::Record)
    return hasseqlength(record) && any(i < 0xff for i in quality(record))
end

"""
    auxdata(record::Record)::BAM.AuxData

Get the auxiliary data of `record`.
"""
function auxdata(record::Record)
    return AuxData(record.data[auxdata_position(record):data_size(record)])
end

function hasauxdata(record::Record)
    return auxdata_position(record) < data_size(record)
end

function Base.getindex(record::Record, tag::AbstractString)
    checkauxtag(tag)
    start = auxdata_position(record)
    stop = data_size(record)
    return getauxvalue(record.data, start, stop, UInt8(tag[1]), UInt8(tag[2]))
end

function Base.haskey(record::Record, tag::AbstractString)
    checkauxtag(tag)
    start = auxdata_position(record)
    stop = data_size(record)
    return findauxtag(record.data, start, stop, UInt8(tag[1]), UInt8(tag[2])) > 0
end

function Base.keys(record::Record)
    return collect(keys(auxdata(record)))
end

function Base.values(record::Record)
    return [record[key] for key in keys(record)]
end


# BioGenerics Methods
# -----------

function BioGenerics.isfilled(record::Record)
    return record.block_size != 0
end

function BioGenerics.seqname(record::Record)
    return tempname(record)
end

function BioGenerics.hasseqname(record::Record)
        return hastempname(record)
end

function BioGenerics.sequence(record::Record)
    return sequence(record)
end

function BioGenerics.hassequence(record::Record)
    return hassequence(record)
end

function BioGenerics.leftposition(record::Record)
    return position(record)
end

function BioGenerics.hasleftposition(record::Record)
    return hasposition(record)
end

function BioGenerics.rightposition(record::Record)
    return rightposition(record)
end

function BioGenerics.hasrightposition(record::Record)
    return hasrightposition(record)
end


# Helper Functions
# ----------------

# Return the size of the `.data` field.
function data_size(record::Record)
    return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
end

function auxdata_position(record::Record)
    seqlen = seqlength(record)
    seqlen === nothing && return 2
    #      name                   null  CIGAR                           DNA              quals 
    return seqname_length(record) + 1 + 4 * n_cigar_op(record, false) + cld(seqlen, 2) + seqlen + 1
end

# Return the length of the read name.
function seqname_length(record::Record)
    return (record.l_read_name % Int) - 1
end
