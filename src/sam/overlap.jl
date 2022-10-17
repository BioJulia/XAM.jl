function GenomicFeatures.eachoverlap(reader::Reader, refname::AbstractString)
    refname = String(refname)
    idx = findfirst(reader.header.metainfo) do m
        SAM.tag(m) == "SQ" && m["SN"] == refname
    end
    refseqlen = parse(Int, reader.header.metainfo[idx]["LN"])
    return GenomicFeatures.eachoverlap(reader, Interval(refname, 1, refseqlen))
end
