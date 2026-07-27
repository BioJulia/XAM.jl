module BGZFStreamsExt

using XAM
using BGZFStreams

import XAM.BAM: write_header

function write_header(stream::BGZFStreams.BGZFStream, header, refseqnames, refseqlens)
    @assert length(refseqnames) == length(refseqlens) "Lengths of refseq names and lengths must match."
    n = 0

    # magic bytes
    n += write(stream, "BAM\1")

    # SAM header
    buf = IOBuffer()
    l = write(SAM.Writer(buf), header)
    n += write(stream, Int32(l))
    n += write(stream, take!(buf))

    # reference sequences
    n += write(stream, Int32(length(refseqnames)))
    for (seqname, seqlen) in zip(refseqnames, refseqlens)
        namelen = length(seqname)
        n += write(stream, Int32(namelen + 1))
        n += write(stream, seqname, '\0')
        n += write(stream, Int32(seqlen))
    end

    return n
end

end