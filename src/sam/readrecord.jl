# Automa.jl generated readrecord! and readmetainfo! functions
# ========================================

# file   = header . body
# header = metainfo*
# body   = record*
const sam_machine_metainfo, sam_machine_record, sam_machine_header, sam_machine_body, sam_machine = (function ()

    isinteractive() && info("compiling SAM")

    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any

    metainfo = let
        tag = re"[A-Z][A-Z]" \ cat("CO")
        tag.actions[:enter] = [:pos1]
        tag.actions[:exit]  = [:metainfo_tag]

        dict = let
            key = re"[A-Za-z][A-Za-z0-9]"
            key.actions[:enter] = [:pos2]
            key.actions[:exit]  = [:metainfo_dict_key]
            val = re"[ -~]+"
            val.actions[:enter] = [:pos2]
            val.actions[:exit]  = [:metainfo_dict_val]
            keyval = cat(key, ':', val)

            cat(keyval, rep(cat('\t', keyval)))
        end
        dict.actions[:enter] = [:pos1]
        dict.actions[:exit]  = [:metainfo_val]

        co = cat("CO")
        co.actions[:enter] = [:pos1]
        co.actions[:exit]  = [:metainfo_tag]

        comment = re"[^\r\n]*"
        comment.actions[:enter] = [:pos1]
        comment.actions[:exit]  = [:metainfo_val]

        cat('@', alt(cat(tag, '\t', dict), cat(co, '\t', comment)))
    end
    metainfo.actions[:enter] = [:mark]
    metainfo.actions[:exit]  = [:metainfo]

    record = let
        qname = re"[!-?A-~]+"
        qname.actions[:enter] = [:pos]
        qname.actions[:exit]  = [:record_qname]

        flag = re"[0-9]+"
        flag.actions[:enter] = [:pos]
        flag.actions[:exit]  = [:record_flag]

        rname = re"\*|[!-()+-<>-~][!-~]*"
        rname.actions[:enter] = [:pos]
        rname.actions[:exit]  = [:record_rname]

        pos = re"[0-9]+"
        pos.actions[:enter] = [:pos]
        pos.actions[:exit]  = [:record_pos]

        mapq = re"[0-9]+"
        mapq.actions[:enter] = [:pos]
        mapq.actions[:exit]  = [:record_mapq]

        cigar = re"\*|([0-9]+[MIDNSHPX=])+"
        cigar.actions[:enter] = [:pos]
        cigar.actions[:exit]  = [:record_cigar]

        rnext = re"\*|=|[!-()+-<>-~][!-~]*"
        rnext.actions[:enter] = [:pos]
        rnext.actions[:exit]  = [:record_rnext]

        pnext = re"[0-9]+"
        pnext.actions[:enter] = [:pos]
        pnext.actions[:exit]  = [:record_pnext]

        tlen = re"[-+]?[0-9]+"
        tlen.actions[:enter] = [:pos]
        tlen.actions[:exit]  = [:record_tlen]

        seq = re"\*|[A-Za-z=.]+"
        seq.actions[:enter] = [:pos]
        seq.actions[:exit]  = [:record_seq]

        qual = re"[!-~]+"
        qual.actions[:enter] = [:pos]
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
        field.actions[:enter] = [:pos]
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
    record.actions[:enter] = [:mark]
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
    body.actions[:exit]  = [:body]

    sam = cat(header, body)

    return map(Automa.compile, (metainfo, record, header, body, sam))
end)()

# write("sam_machine_metainfo.dot", Automa.machine2dot(sam_machine_metainfo))
# run(`dot -Tsvg -o sam_machine_metainfo.svg sam_machine_metainfo.dot`)
#
# write("sam_machine_record.dot", Automa.machine2dot(sam_machine_record))
# run(`dot -Tsvg -o sam_machine_record.svg sam_machine_record.dot`)
#
# write("sam_machine_header.dot", Automa.machine2dot(sam_machine_header))
# run(`dot -Tsvg -o sam_machine_header.svg sam_machine_header.dot`)
#
# write("sam_machine_body.dot", Automa.machine2dot(sam_machine_body))
# run(`dot -Tsvg -o sam_machine_body.svg sam_machine_body.dot`)
#
# write("sam_machine.dot", Automa.machine2dot(sam_machine))
# run(`dot -Tsvg -o sam_machine.svg sam_machine.dot`)

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

const action_metainfo = quote

    let markpos = @markpos()

        appendfrom!(metainfo.data, 1, data, markpos, length(markpos:p-1))

        metainfo.filled = @relpos(markpos):@relpos(p-1)

        found_metainfo = true
    end

end

const sam_actions_metainfo = Dict(
    :mark => :(@mark),
    :pos1  => :(pos1 = @relpos(p)),
    :pos2  => :(pos2 = @relpos(p)),
    :metainfo_tag => :(metainfo.tag = pos1:@relpos(p-1)),
    :metainfo_val => :(metainfo.val = pos1:@relpos(p-1)),
    :metainfo_dict_key => :(push!(metainfo.dictkey, pos2:@relpos(p-1))),
    :metainfo_dict_val => :(push!(metainfo.dictval, pos2:@relpos(p-1))),
    :metainfo => action_metainfo
)

const sam_actions_header = merge(
    sam_actions_metainfo,
    Dict(
        :countline => :(linenum += 1),
        :metainfo => quote
            $(action_metainfo)
            push!(header, metainfo)
            metainfo = MetaInfo()
        end,
        :header => quote

            finish_header = true

            if !eof(stream)
                p -= 1 # cancel look-ahead
            end

            @escape
        end
    )
)

const sam_actions_record = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :record_qname => :(record.qname = pos:@relpos(p-1)),
    :record_flag  => :(record.flag  = pos:@relpos(p-1)),
    :record_rname => :(record.rname = pos:@relpos(p-1)),
    :record_pos   => :(record.pos   = pos:@relpos(p-1)),
    :record_mapq  => :(record.mapq  = pos:@relpos(p-1)),
    :record_cigar => :(record.cigar = pos:@relpos(p-1)),
    :record_rnext => :(record.rnext = pos:@relpos(p-1)),
    :record_pnext => :(record.pnext = pos:@relpos(p-1)),
    :record_tlen  => :(record.tlen  = pos:@relpos(p-1)),
    :record_seq   => :(record.seq   = pos:@relpos(p-1)),
    :record_qual  => :(record.qual  = pos:@relpos(p-1)),
    :record_field => :(push!(record.fields, pos:@relpos(p-1))),
    :record       => quote
        let markpos = @markpos()

            appendfrom!(record.data, 1, data, markpos, length(markpos:p-1))

            record.filled = @relpos(markpos):@relpos(p-1)

            found_record = true
            @escape
        end
    end
)

const sam_actions_body = merge(
    sam_actions_record,
    Dict(
        :countline => :(linenum += 1),
        :body => quote
            finish_body = true
            @escape
        end
    )
)

# const sam_actions = merge(
#     sam_actions_header,
#     sam_actions_body
# )

const sam_context = Automa.CodeGenContext(
    generator = :goto,
    checkbounds = false,
    loopunroll = 0
)

const sam_initcode_metainfo = quote
    pos1 = 0
    pos2 = 0
    found_metainfo = false
end

const sam_initcode_record = quote
    pos = 0
    found_record = false
end

const sam_initcode_header = quote
    $(sam_initcode_metainfo)
    metainfo = MetaInfo()
    finish_header = false
    cs, linenum = state
end

const sam_initcode_body = quote
    $(sam_initcode_record)
    finish_body = false
    cs, linenum = state
end

const sam_loopcode_metainfo = quote

    if cs < 0
        throw(ArgumentError("malformed metainfo at pos $(p)"))
    end

    if found_metainfo
        @goto __return__
    end
end

const sam_returncode_metainfo = quote
    return found_metainfo
end


Automa.Stream.generate_reader(
    :index!,
    sam_machine_metainfo,
    arguments = (:(metainfo::MetaInfo),),
    actions = sam_actions_metainfo,
    context = sam_context,
    initcode = sam_initcode_metainfo,
    loopcode = sam_loopcode_metainfo,
    returncode = sam_returncode_metainfo
) |> eval

const sam_loopcode_header = quote

    if cs < 0
        throw(ArgumentError("malformed metainfo at line $(linenum)"))
    end

    if finish_header
        @goto __return__
    end
end

const sam_returncode_header = quote
    return cs, linenum, finish_header
end

Automa.Stream.generate_reader(
    :readheader!,
    sam_machine_header,
    arguments = (:(header::SAM.Header), :(state::Tuple{Int,Int})),
    actions = sam_actions_header,
    context = sam_context,
    initcode = sam_initcode_header,
    loopcode = sam_loopcode_header,
    returncode = sam_returncode_header
) |> eval


const sam_loopcode_record = quote

    if cs < 0
        throw(ArgumentError("malformed SAM record at position $(p), line $(linenum)"))
    end

    # # if cs != 0
    # #     throw(ArgumentError(string("failed to index ", $(record_type), " ~>", repr(String(data[p:min(p+7,p_end)])))))
    # # end

    if found_record
        @goto __return__
    end

end

const sam_loopcode_body = quote

    $(sam_loopcode_record)

    if finish_body
        @goto __return__
    end
end

const sam_returncode_record = quote
    return found_record
end

Automa.Stream.generate_reader(
    :index!,
    sam_machine_record,
    arguments = (:(record::Record),),
    actions = sam_actions_record,
    context = sam_context,
    initcode = sam_initcode_record,
    loopcode = sam_loopcode_record,
    returncode = sam_returncode_record
) |> eval



const sam_returncode_body = quote
    return cs, linenum, found_record
end

Automa.Stream.generate_reader(
    :readrecord!,
    sam_machine_body,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = sam_actions_body,
    context = sam_context,
    initcode = sam_initcode_body,
    loopcode = sam_loopcode_body,
    returncode = sam_returncode_body
) |> eval
