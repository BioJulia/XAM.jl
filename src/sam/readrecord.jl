@inline function anchor!(stream::BufferedStreams.BufferedInputStream, p, immobilize = true)
    stream.anchor = p
    stream.immobilized = immobilize
    return stream
end

@inline function upanchor!(stream::BufferedStreams.BufferedInputStream)
    @assert stream.anchor != 0 "upanchor! called with no anchor set"
    anchor = stream.anchor
    stream.anchor = 0
    stream.immobilized = false
    return anchor
end

function ensure_margin!(stream::BufferedStreams.BufferedInputStream)
    if stream.position * 20 > length(stream.buffer) * 19
        BufferedStreams.shiftdata!(stream)
    end
    return nothing
end

@inline function resize_and_copy!(dst::Vector{UInt8}, src::Vector{UInt8}, r::UnitRange{Int})
    return resize_and_copy!(dst, 1, src, r)
end

@inline function resize_and_copy!(dst::Vector{UInt8}, dstart::Int, src::Vector{UInt8}, r::UnitRange{Int})
    rlen = length(r)
    if length(dst) != dstart + rlen - 1
        resize!(dst, dstart + rlen - 1)
    end
    copyto!(dst, dstart, src, first(r), rlen)
    return dst
end

function generate_index_function(record_type, machine, init_code, actions; kwargs...)
    kwargs = Dict(kwargs)
    context = Automa.CodeGenContext(
        generator = get(kwargs, :generator, :goto),
        checkbounds = get(kwargs, :checkbounds, false),
        loopunroll = get(kwargs, :loopunroll, 0)
    )
    quote
        function index!(record::$(record_type))
            data = record.data
            p = 1
            p_end = p_eof = sizeof(data)
            initialize!(record)
            $(init_code)
            cs = $(machine.start_state)
            $(Automa.generate_exec_code(context, machine, actions))
            if cs != 0
                throw(ArgumentError(string("failed to index ", $(record_type), " ~>", repr(String(data[p:min(p+7,p_end)])))))
            end
            @assert isfilled(record)
            return record
        end
    end
end

function generate_readheader_function(reader_type, metainfo_type, machine, init_code, actions, finish_code=:())
    quote
        function readheader!(reader::$(reader_type))
            _readheader!(reader, reader.state)
        end

        function _readheader!(reader::$(reader_type), state::State)
            stream = state.stream
            ensure_margin!(stream)
            cs = state.cs
            linenum = state.linenum
            data = stream.buffer
            p = stream.position
            p_end = stream.available
            p_eof = -1
            finish_header = false
            record = $(metainfo_type)()

            $(init_code)

            while true
                $(Automa.generate_exec_code(Automa.CodeGenContext(generator=:table), machine, actions))

                state.cs = cs
                state.finished = cs == 0
                state.linenum = linenum
                stream.position = p

                if cs < 0
                    error("$($(reader_type)) file format error on line ", linenum)
                elseif finish_header
                    $(finish_code)
                    break
                elseif p > p_eof ≥ 0
                    error("incomplete $($(reader_type)) input on line ", linenum)
                else
                    hits_eof = BufferedStreams.fillbuffer!(stream) == 0
                    p = stream.position
                    p_end = stream.available
                    if hits_eof
                        p_eof = p_end
                    end
                end
            end
        end
    end
end

function generate_read_function(reader_type, machine, init_code, actions; kwargs...)
    kwargs = Dict(kwargs)
    context = Automa.CodeGenContext(
        generator=get(kwargs, :generator, :goto),
        checkbounds=get(kwargs, :checkbounds, false),
        loopunroll=get(kwargs, :loopunroll, 0)
    )
    quote
        function Base.read!(reader::$(reader_type), record::eltype($(reader_type)))::eltype($(reader_type))
            return _read!(reader, reader.state, record)
        end

        function _read!(reader::$(reader_type), state::State, record::eltype($(reader_type)))
            stream = state.stream
            ensure_margin!(stream)
            cs = state.cs
            linenum = state.linenum
            data = stream.buffer
            p = stream.position
            p_end = stream.available
            p_eof = -1
            found_record = false
            initialize!(record)

            $(init_code)

            if state.finished
                throw(EOFError())
            end

            while true
                $(Automa.generate_exec_code(context, machine, actions))

                state.cs = cs
                state.finished |= cs == 0
                state.linenum = linenum
                stream.position = p

                if cs < 0
                    error($(reader_type), " file format error on line ", linenum, " ~>", repr(String(data[p:min(p+7,p_end)])))
                elseif found_record
                    break
                elseif cs == 0
                    throw(EOFError())
                elseif p > p_eof ≥ 0
                    error("incomplete $($(reader_type)) input on line ", linenum)
                elseif BufferedStreams.available_bytes(stream) < 64
                    hits_eof = BufferedStreams.fillbuffer!(stream) == 0
                    p = stream.position
                    p_end = stream.available
                    if hits_eof
                        p_eof = p_end
                    end
                end
            end

            @assert isfilled(record)
            return record
        end
    end
end

# Automa.jl generated readrecord! and readmetainfo! functions
# ========================================

# file   = header . body
# header = metainfo*
# body   = record*
const sam_metainfo_machine, sam_record_machine, sam_header_machine, sam_body_machine = (function ()

    isinteractive() && info("compiling SAM")

    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any

    metainfo = let
        tag = re"[A-Z][A-Z]" \ cat("CO")
        tag.actions[:enter] = [:mark1]
        tag.actions[:exit]  = [:metainfo_tag]

        dict = let
            key = re"[A-Za-z][A-Za-z0-9]"
            key.actions[:enter] = [:mark2]
            key.actions[:exit]  = [:metainfo_dict_key]
            val = re"[ -~]+"
            val.actions[:enter] = [:mark2]
            val.actions[:exit]  = [:metainfo_dict_val]
            keyval = cat(key, ':', val)

            cat(keyval, rep(cat('\t', keyval)))
        end
        dict.actions[:enter] = [:mark1]
        dict.actions[:exit]  = [:metainfo_val]

        co = cat("CO")
        co.actions[:enter] = [:mark1]
        co.actions[:exit]  = [:metainfo_tag]

        comment = re"[^\r\n]*"
        comment.actions[:enter] = [:mark1]
        comment.actions[:exit]  = [:metainfo_val]

        cat('@', alt(cat(tag, '\t', dict), cat(co, '\t', comment)))
    end
    metainfo.actions[:enter] = [:anchor]
    metainfo.actions[:exit]  = [:metainfo]

    record = let
        qname = re"[!-?A-~]+"
        qname.actions[:enter] = [:mark]
        qname.actions[:exit]  = [:record_qname]

        flag = re"[0-9]+"
        flag.actions[:enter] = [:mark]
        flag.actions[:exit]  = [:record_flag]

        rname = re"\*|[!-()+-<>-~][!-~]*"
        rname.actions[:enter] = [:mark]
        rname.actions[:exit]  = [:record_rname]

        pos = re"[0-9]+"
        pos.actions[:enter] = [:mark]
        pos.actions[:exit]  = [:record_pos]

        mapq = re"[0-9]+"
        mapq.actions[:enter] = [:mark]
        mapq.actions[:exit]  = [:record_mapq]

        cigar = re"\*|([0-9]+[MIDNSHPX=])+"
        cigar.actions[:enter] = [:mark]
        cigar.actions[:exit]  = [:record_cigar]

        rnext = re"\*|=|[!-()+-<>-~][!-~]*"
        rnext.actions[:enter] = [:mark]
        rnext.actions[:exit]  = [:record_rnext]

        pnext = re"[0-9]+"
        pnext.actions[:enter] = [:mark]
        pnext.actions[:exit]  = [:record_pnext]

        tlen = re"[-+]?[0-9]+"
        tlen.actions[:enter] = [:mark]
        tlen.actions[:exit]  = [:record_tlen]

        seq = re"\*|[A-Za-z=.]+"
        seq.actions[:enter] = [:mark]
        seq.actions[:exit]  = [:record_seq]

        qual = re"[!-~]+"
        qual.actions[:enter] = [:mark]
        qual.actions[:exit]  = [:record_qual]

        field = let
            tag = re"[A-Za-z][A-Za-z0-9]"
            val = alt(
                re"A:[!-~]",
                re"i:[-+]?[0-9]+",
                re"f:[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?",
                re"Z:[ !-~]*",
                re"H:([0-9A-F][0-9A-F])*",
                re"B:[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+")

            cat(tag, ':', val)
        end
        field.actions[:enter] = [:mark]
        field.actions[:exit]  = [:record_field]

        cat(
            qname, '\t',
            flag,  '\t',
            rname, '\t',
            pos,   '\t',
            mapq,  '\t',
            cigar, '\t',
            rnext, '\t',
            pnext, '\t',
            tlen,  '\t',
            seq,   '\t',
            qual,
            rep(cat('\t', field)))
    end
    record.actions[:enter] = [:anchor]
    record.actions[:exit]  = [:record]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(re"\r?", lf)
    end

    header′ = rep(cat(metainfo, newline))
    header′.actions[:exit] = [:header]
    header = cat(header′, opt(any() \ cat('@')))  # look ahead

    body = rep(cat(record, newline))

    return map(Automa.compile, (metainfo, record, header, body))
end)()

const sam_metainfo_actions = Dict(
    :metainfo_tag => :(record.tag = (mark1:p-1) .- offset),
    :metainfo_val => :(record.val = (mark1:p-1) .- offset),
    :metainfo_dict_key => :(push!(record.dictkey, (mark2:p-1) .- offset)),
    :metainfo_dict_val => :(push!(record.dictval, (mark2:p-1) .- offset)),
    :metainfo => quote
        resize_and_copy!(record.data, data, offset+1:p-1)
        record.filled = (offset+1:p-1) .- offset
    end,
    :anchor => :(),
    :mark1  => :(mark1 = p),
    :mark2  => :(mark2 = p)
)

generate_index_function(
    MetaInfo,
    sam_metainfo_machine,
    :(mark1 = mark2 = offset = 0),
    sam_metainfo_actions
) |> eval

generate_readheader_function(
    Reader,
    MetaInfo,
    sam_header_machine,
    :(mark1 = mark2 = offset = 0),
    merge(sam_metainfo_actions, Dict(
        :metainfo => quote
            resize_and_copy!(record.data, data, upanchor!(stream):p-1)
            record.filled = (offset+1:p-1) .- offset
            @assert isfilled(record)
            push!(reader.header.metainfo, record)
            ensure_margin!(stream)
            record = MetaInfo()
        end,
        :header => :(finish_header = true; @escape),
        :countline => :(linenum += 1),
        :anchor => :(anchor!(stream, p); offset = p - 1))),
    quote
        if !eof(stream)
            stream.position -= 1  # cancel look-ahead
        end
    end
) |> eval

const sam_record_actions = Dict(
    :record_qname => :(record.qname = (mark:p-1) .- offset),
    :record_flag  => :(record.flag  = (mark:p-1) .- offset),
    :record_rname => :(record.rname = (mark:p-1) .- offset),
    :record_pos   => :(record.pos   = (mark:p-1) .- offset),
    :record_mapq  => :(record.mapq  = (mark:p-1) .- offset),
    :record_cigar => :(record.cigar = (mark:p-1) .- offset),
    :record_rnext => :(record.rnext = (mark:p-1) .- offset),
    :record_pnext => :(record.pnext = (mark:p-1) .- offset),
    :record_tlen  => :(record.tlen  = (mark:p-1) .- offset),
    :record_seq   => :(record.seq   = (mark:p-1) .- offset),
    :record_qual  => :(record.qual  = (mark:p-1) .- offset),
    :record_field => :(push!(record.fields, (mark:p-1) .- offset)),
    :record       => quote
        resize_and_copy!(record.data, data, 1:p-1)
        record.filled = (offset+1:p-1) .- offset
    end,
    :anchor       => :(),
    :mark         => :(mark = p)
)

generate_index_function(
    Record,
    sam_record_machine,
    :(mark = offset = 0),
    sam_record_actions
) |> eval

generate_read_function(
    Reader,
    sam_body_machine,
    :(mark = offset = 0),
    merge(sam_record_actions, Dict(
        :record    => quote
            resize_and_copy!(record.data, data, upanchor!(stream):p-1)
            record.filled = (offset+1:p-1) .- offset
            found_record = true
            @escape
        end,
        :countline => :(linenum += 1),
        :anchor    => :(anchor!(stream, p); offset = p - 1))
    )
) |> eval
