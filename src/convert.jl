module Conversion

using ..SAM
using ..BAM

using BioSequences, BioAlignments

global  const op_to_code = Dict(
    'M' => 0x00,
    'I' => 0x01,
    'D' => 0x02,
    'N' => 0x03,
    'S' => 0x04,
    'H' => 0x05,
    'P' => 0x06,
    '=' => 0x07,
    'X' => 0x08,
    'B' => 0x09,
    '0' => 0x0a
)

function parse_cigar(record::SAM.Record)
    
    cigar_ops = ('M','I','D','N','S','H','P','=','X')
    out = NamedTuple{(:L, :op), Tuple{Int64, Char}}[]

    cigar = SAM.cigar(record)
    start = 1

    for (i,c) in enumerate(cigar)
        if c ∈ cigar_ops
            L = parse(Int, cigar[start:i-1])
            op = c
            push!(out, (L = L, op = op))
            start = i+1
        end
    end
    out
end

function Base.convert(::Type{BAM.Record}, record::SAM.Record; header = nothing)

    bam_record = BAM.Record()
    set_fixed_length_fields!(bam_record, record, header)
    set_variable_length_fields!(bam_record, record)

    bam_record
end

function refname_to_refid(refname::String, header)
    isnothing(header) && return -1
    refname == "*" && return -1
    
    refs = findall(header, "SQ")
    for (i,ref) in enumerate(refs)
        ref["SN"] == refname && return i-1
    end
    error("Couldn't find reference $(refname) in header.")
end

"""
    This field is set as '*' when the information is unavailable, and set as '=' if RNEXT is identical RNAME
"""
function get_nextrefname(record::SAM.Record)
    !SAM.hasnextrefname(record) && return "*"
    SAM.nextrefname(record) == "=" ? SAM.refname(record) : SAM.nextrefname(record)
end

function get_refname(record::SAM.Record)
    !SAM.hasrefname(record) && return "*"
    SAM.refname(record)
end

function set_fixed_length_fields!(bam_record::BAM.Record, record::SAM.Record, header)

    cigar_ops = ('M','I','D','N','S','H','P','=','X')
    n_cigar_op = sum(x ∈ cigar_ops for x in SAM.cigar(record); init = 0)

    if n_cigar_op > typemax(UInt16) 
        error("Cigar contains more than $(typemax(UInt16)) operations. This is currently unsupported. See section 4.2.2 of SAM/BAM specifications.")
    end

    bin = reg2bin(SAM.position(record)-1, SAM.alignlength(record))
    
    bam_record.refid = refname_to_refid(get_refname(record), header)
    bam_record.pos = SAM.position(record)-1
    bam_record.l_read_name = SAM.tempname(record) |> x -> length(x)+1 # +1 for NULL
    bam_record.mapq = SAM.mappingquality(record)
    bam_record.bin = bin
    bam_record.n_cigar_op = n_cigar_op
    bam_record.flags = SAM.flags(record)
    bam_record.l_seq = SAM.hassequence(record) ? SAM.seqlength(record) : 0
    bam_record.next_refid = refname_to_refid(get_nextrefname(record), header)
    bam_record.next_pos = SAM.nextposition(record)-1
    bam_record.tlen = SAM.templength(record)
    bam_record
end

"""
While all single (i.e., non-array) integer types are stored in SAM as 'i', in BAM any of 'cCsSiI' 
may be used together with the correspondingly-sized binary integer value, chosen according to the field value's magnitude
"""
global const type_to_aux_int_fields = Dict(
    Int8 => UInt8('c'),
    UInt8 => UInt8('C'),
    Int16 => UInt8('s'),
    UInt16 => UInt8('S'),
    Int32 => UInt8('i'),
    UInt32 => UInt8('I'),
)

function find_smallest_int_type(x)
    candidates = x < 0 ? (Int8, Int16, Int32) : (UInt8, UInt16, UInt32)
    idx = findfirst(T -> typemin(T) <= x <= typemax(T) , candidates)
    isnothing(idx) && error("Coulnd't find Int type for value $(x)")
    T = candidates[idx]

    T, type_to_aux_int_fields[T]
end

function encode_aux_field(record::SAM.Record, field::UnitRange{Int64})
    typ = record.data[first(field)+3]
    key = record.data[first(field):first(field)+1]

    if typ == UInt8('i')
        val = record.data[first(field)+5:last(field)] |> String
        val = parse(Int32, val)
        T, typ = find_smallest_int_type(val)
        val = reinterpret(UInt8, [T(val)])

        return vcat(key, typ, val)
    end
    if typ == UInt8('Z')
        val = record.data[first(field)+5:last(field)]
        return vcat(key, typ, val, 0x00)
    end
    if typ == UInt8('A')
        val = record.data[first(field)+5:last(field)]
        return vcat(key, typ, val)
    end

    error("Unsupported tag type $(Char(typ))")
end

function encode_cigar(record::SAM.Record)

    ops = parse_cigar(record)
    cigar = Vector{UInt8}(undef, length(ops)*4)
    k = 1
    for c in ops
        cigar_32 = [UInt32(0) + op_to_code[c.op] + (UInt32(c.L) << 4)]
        for u in reinterpret(UInt8, cigar_32)
            cigar[k] = u
            k+=1
        end
    end
    cigar
end

function encode_sequence(record::SAM.Record)
    !SAM.hassequence(record) && return UInt8[]

    Nsam = SAM.seqlength(record)
    Nbam = iseven(Nsam) ? div(Nsam,2) : div(Nsam+1,2)

    # Specs : "When l seq is odd the bottom 4 bits of the last byte are undefined, but we recommend writing these as zero."
    bamseq = fill(0x00, Nbam)
    samseq = SAM.sequence(record)

    # SAM
    # 1 2 3 4 5
    # A T G G T
    
    # BAM
    # 1  2  3
    # AT GG T0

    for i in 1:Nbam
        k = 2i -1
        val = BioSequences.encode(DNAAlphabet{4}(), samseq[k]) |> UInt8
        bamseq[i] += val <<4
        
        k = 2i
        k > Nsam && break
        val = BioSequences.encode(DNAAlphabet{4}(), samseq[k]) |> UInt8
        bamseq[i] += val
    end

    bamseq
end

function encode_tempname(record::SAM.Record)
    tempname = SAM.tempname(record)
    vcat(unsafe_wrap(Vector{UInt8}, tempname), 0x00) #add NULL terminator
end

""" See section 5.3 of specifications."""
function reg2bin(start::Int, stop::Int)
    stop -= 1
    (start >> 14) == (stop >> 14) && return ((1 << 15) - 1) ÷ 7 + (start >> 14)
    (start >> 17) == (stop >> 17) && return ((1 << 12) - 1) ÷ 7 + (start >> 17)
    (start >> 20) == (stop >> 20) && return ((1 << 9)  - 1) ÷ 7 + (start >> 20)
    (start >> 23) == (stop >> 23) && return ((1 << 6)  - 1) ÷ 7 + (start >> 23)
    (start >> 26) == (stop >> 26) && return ((1 << 3)  - 1) ÷ 7 + (start >> 26)
    return 0
end

function set_variable_length_fields!(bam_record::BAM.Record, record::SAM.Record)

    tempname = encode_tempname(record)
    seq = encode_sequence(record)
    cigar = encode_cigar(record)
    aux_data = reduce(vcat, encode_aux_field(record, field) for field in record.fields; init = UInt8[])
    # if quality is missing, fill with good quality
    quality = SAM.hasquality(record) ? SAM.quality(record) : fill(0xff, SAM.hassequence(record) ? SAM.seqlength(record) : 0) 

    bam_record.data = vcat(
        tempname,
        cigar,
        seq,
        quality,
        aux_data
    )
    bam_record.block_size = length(bam_record.data) + BAM.FIXED_FIELDS_BYTES - sizeof(bam_record.block_size)
    bam_record
end

end# module Conversion