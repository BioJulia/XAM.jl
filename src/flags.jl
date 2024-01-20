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
        (:PAIRED,        UInt16(0x001), "the read is paired in sequencing, no matter whether it is mapped in a pair"),
        (:PROPER_PAIR,   UInt16(0x002), "the read is mapped in a proper pair"                                       ),
        (:UNMAP,         UInt16(0x004), "the read itself is unmapped; conflictive with FLAG_PROPER_PAIR"            ),
        (:MUNMAP,        UInt16(0x008), "the mate is unmapped"                                                      ),
        (:REVERSE,       UInt16(0x010), "the read is mapped to the reverse strand"                                  ),
        (:MREVERSE,      UInt16(0x020), "the mate is mapped to the reverse strand"                                  ),
        (:READ1,         UInt16(0x040), "this is read1"                                                             ),
        (:READ2,         UInt16(0x080), "this is read2"                                                             ),
        (:SECONDARY,     UInt16(0x100), "not primary alignment"                                                     ),
        (:QCFAIL,        UInt16(0x200), "QC failure"                                                                ),
        (:DUP,           UInt16(0x400), "optical or PCR duplicate"                                                  ),
        (:SUPPLEMENTARY, UInt16(0x800), "supplementary alignment"                                                   ),]
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

Query whether the `record`'s template has multiple segments in sequencing.
"""
function ispaired(record::XAMRecord)::Bool
    return flags(record) & FLAG_PAIRED == FLAG_PAIRED
end

"""
    isproperpair(record::XAMRecord)::Bool

Query whether each segment of the `record`'s template properly aligned according to the aligner.
"""
function isproperpair(record::XAMRecord)::Bool
    return flags(record) & PROPER_PAIR == PROPER_PAIR
end

"""
    isunmapped(record::XAMRecord)::Bool

Query whether the `record` is unmapped.
"""
function isunmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_UNMAP == FLAG_UNMAP
end

"""
    ismapped(record::XAMRecord)::Bool

Query whether the `record` is mapped.
"""
function ismapped(record::XAMRecord)::Bool
    # return flags(record) & FLAG_UNMAP == 0
    return isfilled(record) && (flags(record) & FLAG_UNMAP == 0)
end

"""
    ismateunmapped(record::XAMRecord)::Bool

Query whether the `record`'s mate is unmapped.
"""
function ismateunmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_MUNMAP == FLAG_MUNMAP
end

"""
    isnextmapped(record::XAMRecord)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::XAMRecord)::Bool
    return flags(record) & FLAG_MUNMAP == 0
end

"""
    isreverse(record::XAMRecord)::Bool

Query whether the `record` is mapped to the reverse strand.
"""
function isreverse(record::XAMRecord)::Bool
    return flags(record) & FLAG_REVERSE == FLAG_REVERSE
end

"""
    isforward(record::XAMRecord)::Bool

Query whether the `record` is mapped to the forward strand.
"""
function isforward(record::XAMRecord)::Bool
    return flags(record) & FLAG_REVERSE == 0
end

"""
    ispositivestrand(record::XAMRecord)::Bool

Query whether `record` is aligned to the positive strand.
"""
function ispositivestrand(record::XAMRecord)::Bool
    return isforward(record)
end

"""
    ispositivestrand(record::XAMRecord)::Bool

Query whether `record` is aligned to the negative strand.
"""
function isnegativestrand(record::XAMRecord)::Bool
    return isreverse(record)
end

"""
    ismatereverse(record::XAMRecord)::Bool

Query whether the `record`'s mate is mapped to the reverse strand.
"""
function ismatereverse(record::XAMRecord)::Bool
    return flags(record) & FLAG_MREVERSE == FLAG_MREVERSE
end

"""
    isread1(record::XAMRecord)::Bool

Query whether the `record` is read1.
"""
function isread1(record::XAMRecord)::Bool
    return flags(record) & FLAG_READ1 == FLAG_READ1
end

"""
    isread2(record::XAMRecord)::Bool

Query whether the `record` is read2.
"""
function isread2(record::XAMRecord)::Bool
    return flags(record) & FLAG_READ2 == FLAG_READ2
end

"""
    issecondaryalignment(record::XAMRecord)::Bool

Query whether the `record` is a secondary alignment.
"""
function issecondaryalignment(record::XAMRecord)::Bool
    return flags(record) & FLAG_SECONDARY == FLAG_SECONDARY
end

"""
    isprimaryalignment(record::XAMRecord)::Bool

Query whether the `record` is the primary alignment.
"""
function isprimaryalignment(record::XAMRecord)::Bool
    return flags(record) & FLAG_SECONDARY == 0
end

"""
    isqcfail(record::XAMRecord)::Bool

Query whether the `record` did not pass filters, such as platform/vendor quality controls.
"""
function isqcfail(record::XAMRecord)::Bool
    return flags(record) & FLAG_QCFAIL == FLAG_QCFAIL
end

"""
    isduplicate(record::XAMRecord)::Bool

Query whether the `record` is a PCR or optical duplicate.
"""
function isduplicate(record::XAMRecord)::Bool
    return flags(record) & FLAG_DUP == FLAG_DUP
end

"""
    issupplementaryalignment(record::XAMRecord)::Bool

Query whether the `record` is a supplementary alignment.
"""
function issupplementaryalignment(record::XAMRecord)::Bool
    return flags(record) & FLAG_SUPPLEMENTARY == FLAG_SUPPLEMENTARY
end

"""
    isprimary(record::XAMRecord)::Bool

Query whether `record` is a primary line of the read.

This is equivalent to `flags(record) & 0x900 == 0`.
"""
function isprimary(record::XAMRecord)::Bool
    return flags(record) & 0x900 == 0
end
