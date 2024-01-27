# Flags
# =========
#

"""
    flags(record::XAMRecord})::UInt16

Get the bitwise flags of `record`.
The returned value is a `UInt16` of each flag being OR'd together.

The possible flags are:

    0x0001 template having multiple segments in sequencing
    0x0002 each segment properly aligned according to the aligner
    0x0004 segment unmapped
    0x0008 next segment in the template unmapped
    0x0010 SEQ being reverse complemented
    0x0020 SEQ of the next segment in the template being reverse complemented
    0x0040 the first segment in the template
    0x0080 the last segment in the template
    0x0100 secondary alignment
    0x0200 not passing filters, such as platform/vendor quality controls
    0x0400 PCR or optical duplicate
    0x0800 supplementary alignment
"""

function flags end




# Bitwise flag (or FLAG).
for (name, bits, doc) in [
        (:PAIRED,          UInt16(0x001), "the segment is paired with other segments"),
        (:PROPER_PAIR,     UInt16(0x002), "the segment is in a template where all segments are properly aligned according to the aligner"),
        (:UNMAPPED,        UInt16(0x004), "the segment itself is unmapped; conflictive with FLAG_PROPER_PAIR"),
        (:NEXT_UNMAPPED,   UInt16(0x008), "the next segment in the template is unmapped"),
        (:REVERSE,         UInt16(0x010), "the *SEQ*uence is reverse complemented"),
        (:NEXT_REVERSE,    UInt16(0x020), "the *SEQ*uence of the next segment in the template is reverse complemented"                                  ),
        (:FIRST_SEGMENT,   UInt16(0x040), "the segment is the first in the template"),
        (:LAST_SEGMENT,    UInt16(0x080), "the segment is last in the template"),
        (:SECONDARY,       UInt16(0x100), "not primary alignment"),
        (:QCFAIL,          UInt16(0x200), "QC failure"),
        (:DUPLICATE,       UInt16(0x400), "optical or PCR duplicate"),
        (:SUPPLEMENTARY,   UInt16(0x800), "supplementary alignment"),
    ]
    @assert bits isa UInt16 "The bits must be of type UInt16."
    sym = Symbol("FLAG_", name)
    docstring = """    $sym
    SAM/BAM flags: $doc

    See also: [`flags`](@ref)
    """
    @eval begin
        @doc $(docstring) const $(sym) = $(bits)
    end
end

"""
    ispaired(record::XAMRecord)::Bool

Query whether the segment is in a template having multiple segments in sequencing.
"""
function ispaired(record::XAMRecord)::Bool
    return flags(record) & FLAG_PAIRED == FLAG_PAIRED
end

"""
    isproperpair(record::XAMRecord)::Bool

Query whether the segment is in a template where all segments are properly aligned according to the aligner.
"""
function isproperpair(record::XAMRecord)::Bool
    return flags(record) & PROPER_PAIR == PROPER_PAIR
end

"""
    isunmapped(record::XAMRecord)::Bool

Query whether the segment is unmapped.
"""
function isunmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_UNMAPPED == FLAG_UNMAPPED
end

"""
    ismapped(record::XAMRecord)::Bool

Query whether the segment is mapped.
"""
function ismapped(record::XAMRecord)::Bool
    # return flags(record) & FLAG_UNMAPPED == 0
    return isfilled(record) && (flags(record) & FLAG_UNMAPPED == 0)
end

"""
    isnextunmapped(record::XAMRecord)::Bool

Query whether the next segment in the template is unmapped.
"""
function isnextunmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_NEXT_UNMAPPED == FLAG_NEXT_UNMAPPED
end

"""
    isnextmapped(record::XAMRecord)::Bool

Query whether the next segment in the template is mapped.
"""
function isnextmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_NEXT_UNMAPPED == 0
end

"""
    isreverse(record::XAMRecord)::Bool

Query whether the `record.SEQ`uence is reverse complemented.
"""
function isreversecomplemented(record::XAMRecord)::Bool
    return flags(record) & FLAG_REVERSE == FLAG_REVERSE
end

"""
    isforward(record::XAMRecord)::Bool

Query whether the `record.SEQ`uence is mapped to the forward strand.
"""
function isforwardstrand(record::XAMRecord)::Bool
    # return flags(record) & FLAG_REVERSE == 0
    return !isreversecomplemented(record) # Note: this is an interpretation of FLAG_REVERSE.
end

"""
    ispositivestrand(record::XAMRecord)::Bool

Query whether the `record.SEQ`uence is aligned to the positive strand.
"""
function ispositivestrand(record::XAMRecord)::Bool
    return isforwardstrand(record)
end

"""
    isreversestrand(record::XAMRecord)::Bool

Query whether the `record.SEQ`uence is aligned to the reverse strand.
"""
function isreversestrand(record::XAMRecord)::Bool
    return isreversecomplemented(record) # Note: this is an interpretation of FLAG_REVERSE.
end

"""
    ispositivestrand(record::XAMRecord)::Bool

Query whether the `record.SEQ`uence is aligned to the negative strand.
"""
function isnegativestrand(record::XAMRecord)::Bool
    return isreversestrand(record)
end

"""
    isnextreversecomplemented(record::XAMRecord)::Bool

Query whether the next segment in the template is reverse complemented.
"""
function isnextreversecomplemented(record::XAMRecord)::Bool
    return flags(record) & FLAG_NEXT_REVERSE == FLAG_NEXT_REVERSE
end

"""
    isfirstsegment(record::XAMRecord)::Bool

Query whether the segemnt is first in the template.
"""
function isfirstsegment(record::XAMRecord)::Bool
    return flags(record) & FLAG_FIRST_SEGMENT == FLAG_FIRST_SEGMENT
end

"""
    isread1(record::XAMRecord)::Bool

From a paired-end sequencing point of view, query whether the read is read1.
"""
function isread1(record::XAMRecord)::Bool
    return isfirstsegment(record)
end

"""
    islastsegment(record::XAMRecord)::Bool

Query whether the segemnt is last in the template.
"""
function islastsegment(record::XAMRecord)::Bool
    return flags(record) & FLAG_LAST_SEGMENT == FLAG_LAST_SEGMENT
end

"""
    isread2(record::XAMRecord)::Bool

From a paired-end sequencing point of view, query whether the read is read2.
"""
function isread2(record::XAMRecord)::Bool
    return islastsegment(record)
end

"""
    issecondaryalignment(record::XAMRecord)::Bool

Query whether the read is considered to be the secondary alignment.
"""
function issecondaryalignment(record::XAMRecord)::Bool
    return flags(record) & FLAG_SECONDARY == FLAG_SECONDARY
end


"""
    isqcfail(record::XAMRecord)::Bool

Query whether the read failed filters, such as platform/vendor quality controls.
"""
function isqcfail(record::XAMRecord)::Bool
    return flags(record) & FLAG_QCFAIL == FLAG_QCFAIL
end

"""
    isduplicate(record::XAMRecord)::Bool

Query whether the read is a PCR or optical duplicate.
"""
function isduplicate(record::XAMRecord)::Bool
    return flags(record) & FLAG_DUPLICATE == FLAG_DUPLICATE
end

"""
    issupplementaryalignment(record::XAMRecord)::Bool

Query whether the read alignment is considered to be a supplementary alignment.
"""
function issupplementaryalignment(record::XAMRecord)::Bool
    return flags(record) & FLAG_SUPPLEMENTARY == FLAG_SUPPLEMENTARY
end

"""
    isprimaryalignment(record::XAMRecord)::Bool

Query whether the read alignment is considered to be the primary alignment.
This is primary line of the read and is equivalent to `flags(record) & 0x900 == 0`.
"""
function isprimaryalignment(record::XAMRecord)::Bool
    # return !issecondaryalignment(record) && !issupplementaryalignment(record)
    return flags(record) & 0x900 == 0
end
