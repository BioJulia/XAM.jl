# Automa.jl generated readrecord! and readmetainfo! functions
# ========================================

import Automa: @re_str, rep, onexit!, onenter!, CodeGenContext, generate_reader, compile

# file   = header . body
# header = metainfo*
# body   = record*
const sam_machine_metainfo, sam_machine_record, sam_machine_header, sam_machine_body, sam_machine = let
    metainfo = let
        tag = onexit!(onenter!(re"[A-Z][A-Z]" \ re"CO", :pos1), :metainfo_tag)

        dict = let
            key = onexit!(onenter!(re"[A-Za-z][A-Za-z0-9]", :pos2), :metainfo_dict_key)
            val = onexit!(onenter!(re"[ -~]+", :pos2), :metainfo_dict_val)
            keyval = key * ':' * val
            keyval * rep('\t' * keyval)
        end
        onexit!(onenter!(dict, :pos1), :metainfo_val)

        co = onexit!(onenter!(re"CO", :pos1), :metainfo_tag)

        comment = onexit!(onenter!(re"[^\r\n]*", :pos1), :metainfo_val) # Note: Only single line comments are allowed.

        '@' * ((tag * '\t' * dict) | (co * '\t' * comment))
    end
    onexit!(onenter!(metainfo, :mark), :metainfo)

    record = let
        qname = onexit!(onenter!(re"[!-?A-~]+", :pos), :record_qname)
        flag = onexit!(onenter!(re"[0-9]+", :pos), :record_flag)
        rname = onexit!(onenter!(re"\*|[!-()+-<>-~][!-~]*", :pos), :record_rname)
        pos = onexit!(onenter!(re"[0-9]+", :pos), :record_pos)
        mapq = onexit!(onenter!(re"[0-9]+", :pos), :record_mapq)
        cigar = onexit!(onenter!(re"\*|([0-9]+[MIDNSHPX=])+", :pos), :record_cigar)
        rnext = onexit!(onenter!(re"\*|=|[!-()+-<>-~][!-~]*", :pos), :record_rnext)
        pnext = onexit!(onenter!(re"[0-9]+", :pos), :record_pnext)
        tlen = onexit!(onenter!(re"[-+]?[0-9]+", :pos), :record_tlen)
        seq = onexit!(onenter!(re"\*|[A-Za-z=.]+", :pos), :record_seq)
        qual = onexit!(onenter!(re"[!-~]+", :pos), :record_qual)

        field = let
            tag = re"[A-Za-z][A-Za-z0-9]"
            val = (
                re"A:[!-~]" |
                re"i:[-+]?[0-9]+" |
                re"f:[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?" |
                re"Z:[ !-~]*" |
                re"H:([0-9A-F][0-9A-F])*" |
                re"B:[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+"
            )
            tag * ':' * val
        end
        onexit!(onenter!(field, :pos), :record_field)

        qname * '\t' *
        flag  * '\t' *
        rname * '\t' *
        pos   * '\t' *
        mapq  * '\t' *
        cigar * '\t' *
        rnext * '\t' *
        pnext * '\t' *
        tlen  * '\t' *
        seq   * '\t' *
        qual *
        rep('\t' * field)
    end
    onexit!(onenter!(record, :mark), :record)
    newline = "\r?" * onenter!(re"\n", :countline)
    header = onexit!(rep(metainfo * newline), :header)
    body = onexit!(rep(record * newline), :body)
    sam = header * body

    map(compile, (metainfo, record, header, body, sam))
end

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

const sam_actions_metainfo = Dict(
    :mark => :(@mark),
    :pos1  => :(pos1 = @relpos(p)),
    :pos2  => :(pos2 = @relpos(p)),
    :metainfo_tag => :(metainfo.tag = pos1:@relpos(p-1)),
    :metainfo_val => :(metainfo.val = pos1:@relpos(p-1)),
    :metainfo_dict_key => :(push!(metainfo.dictkey, pos2:@relpos(p-1))),
    :metainfo_dict_val => :(push!(metainfo.dictval, pos2:@relpos(p-1))),
    :metainfo => quote
        appendfrom!(metainfo.data, 1, data, @markpos, p-@markpos)
        metainfo.filled = 1:(p-@markpos)
    end
)

const sam_actions_header = merge(
    sam_actions_metainfo,
    Dict(
        :countline => :(linenum += 1),
        :metainfo => quote
            $(sam_actions_metainfo[:metainfo])
            push!(header, metainfo)
            metainfo = MetaInfo()
        end,
        :header => :(@escape)
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
        appendfrom!(record.data, 1, data, @markpos, p-@markpos)
        record.filled = 1:(p-@markpos)
    end
)

const sam_actions_body = merge(
    sam_actions_record,
    Dict(
        :countline => :(linenum += 1),
        :record => quote
            found_record = true
            $(sam_actions_record[:record])
            @escape
        end,
        :body => :(@escape)
    )
)

const sam_context = CodeGenContext(generator = :goto)

const sam_initcode_metainfo = quote
    pos1 = 0
    pos2 = 0
end

const sam_initcode_header = quote
    $(sam_initcode_metainfo)
    metainfo = MetaInfo()
    cs, linenum = state
end

const sam_initcode_record = quote
    pos = 0
end

const sam_initcode_body = quote
    $(sam_initcode_record)
    found_record = false
    cs, linenum = state
end

generate_reader(
    :index!,
    sam_machine_metainfo,
    arguments = (:(metainfo::MetaInfo),),
    actions = sam_actions_metainfo,
    context = sam_context,
    initcode = sam_initcode_metainfo,
) |> eval

const sam_returncode_header = quote
    return cs, linenum
end

generate_reader(
    :readheader!,
    sam_machine_header,
    arguments = (:(header::SAM.Header), :(state::Tuple{Int,Int})),
    actions = sam_actions_header,
    context = sam_context,
    initcode = sam_initcode_header,
    returncode = sam_returncode_header,
    errorcode = quote
        # We expect the SAM header machine to error, as it finds the first non-header byte.
        # This happens at state 1 (hence error state -1), and before reaching EOF.
        if cs == -1 && !(is_eof && p < p_end)
            @goto __return__
        else
            error("Expected input byte after SAM header.")
        end
    end
) |> eval

const sam_loopcode_body = quote
    if found_record
        @goto __return__
    end
end

generate_reader(
    :index!,
    sam_machine_record,
    arguments = (:(record::Record),),
    actions = sam_actions_record,
    context = sam_context,
    initcode = sam_initcode_record,
) |> eval


const sam_returncode_body = quote
    return cs, linenum, found_record
end

generate_reader(
    :readrecord!,
    sam_machine_body,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = sam_actions_body,
    context = sam_context,
    initcode = sam_initcode_body,
    loopcode = sam_loopcode_body,
    returncode = sam_returncode_body
) |> eval
