# BAM Overlap
# ===========

struct OverlapIterator{T}
    reader::Reader{T}
    refname::String
    interval::UnitRange{Int}
end

function Base.IteratorSize(::Type{OverlapIterator{T}}) where T
    return Base.SizeUnknown()
end

function Base.eltype(::Type{OverlapIterator{T}}) where T
    return Record
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::Interval)
    return GenomicFeatures.eachoverlap(reader, GenomicFeatures.seqname(interval), GenomicFeatures.leftposition(interval):GenomicFeatures.rightposition(interval))
end

function GenomicFeatures.eachoverlap(reader::Reader, interval)
    return GenomicFeatures.eachoverlap(reader, convert(Interval, interval))
end

function GenomicFeatures.eachoverlap(reader::Reader, refname::AbstractString, interval::UnitRange)
    return OverlapIterator(reader, String(refname), interval)
end


# Iterator
# --------

mutable struct OverlapIteratorState
    # reference index
    refindex::Int

    # possibly overlapping chunks
    chunks::Vector{Indexes.Chunk}

    # current chunk index
    chunkid::Int

    # pre-allocated record
    record::Record
end

function Base.iterate(iter::OverlapIterator)
    refindex = findfirst(isequal(iter.refname), iter.reader.refseqnames)
    if refindex === nothing
        throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    end

    @assert iter.reader.index !== nothing "Reader index cannot be nothing."

    chunks = Indexes.overlapchunks(iter.reader.index.index, refindex, iter.interval)
    if isempty(chunks)
        return nothing
    end
    state = OverlapIteratorState(refindex, chunks, 1, Record())
    seek(iter.reader, state.chunks[state.chunkid].start)
    return iterate(iter, state)
end

function Base.iterate(iter::OverlapIterator, state)
    while state.chunkid ≤ lastindex(state.chunks)
        chunk = state.chunks[state.chunkid]
        while BGZFStreams.virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, state.record)
            c = compare_intervals(state.record, (state.refindex, iter.interval))
            if c == 0  # overlapping
                return copy(state.record), state
            end
            if c > 0
                # no more overlapping records in this chunk since records are sorted
                break
            end
        end
        state.chunkid += 1
        if state.chunkid ≤ lastindex(state.chunks)
            seek(iter.reader, state.chunks[state.chunkid].start)
        end
    end
    # no more overlapping records
    return nothing
end

function compare_intervals(record::Record, interval::Tuple{Int,UnitRange{Int}})
    rid = refid(record)

    if rid < interval[1] || (rid == interval[1] && rightposition(record) < first(interval[2]))
        # strictly left
        return -1
    end

    if rid > interval[1] || (rid == interval[1] && leftposition(record) > last(interval[2]))
        # strictly right
        return +1
    end

    # overlapping
    return 0

end
