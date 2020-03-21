using Test

using BioGenerics
using FormatSpecimens
using GenomicFeatures
using XAM

import BioAlignments: Alignment, AlignmentAnchor, OP_START, OP_MATCH, OP_DELETE
import BGZFStreams: BGZFStream
import BioGenerics.Exceptions: MissingFieldException
import BioSequences: @dna_str, @aa_str


# Generate a random range within `range`.
function randrange(range)
    x = rand(range)
    y = rand(range)
    if x < y
        return x:y
    else
        return y:x
    end
end

@testset "SAM" begin
    samdir = path_of_format("SAM")

    @testset "MetaInfo" begin
        metainfo = SAM.MetaInfo()
        @test !isfilled(metainfo)
        @test occursin("not filled", repr(metainfo))

        metainfo = SAM.MetaInfo("CO", "some comment (parens)")
        @test isfilled(metainfo)
        @test string(metainfo) == "@CO\tsome comment (parens)"
        @test occursin("CO", repr(metainfo))
        @test SAM.tag(metainfo) == "CO"
        @test SAM.value(metainfo) == "some comment (parens)"
        @test_throws ArgumentError keys(metainfo)
        @test_throws ArgumentError values(metainfo)

        metainfo = SAM.MetaInfo("HD", ["VN" => "1.0", "SO" => "coordinate"])
        @test isfilled(metainfo)
        @test string(metainfo) == "@HD\tVN:1.0\tSO:coordinate"
        @test occursin("HD", repr(metainfo))
        @test SAM.tag(metainfo) == "HD"
        @test SAM.value(metainfo) == "VN:1.0\tSO:coordinate"
        @test keys(metainfo) == ["VN", "SO"]
        @test values(metainfo) == ["1.0", "coordinate"]
        @test SAM.keyvalues(metainfo) == ["VN" => "1.0", "SO" => "coordinate"]
        @test haskey(metainfo, "VN")
        @test haskey(metainfo, "SO")
        @test !haskey(metainfo, "GO")
        @test metainfo["VN"] == "1.0"
        @test metainfo["SO"] == "coordinate"
        @test_throws KeyError metainfo["GO"]
    end

    @testset "Header" begin
        header = SAM.Header()
        @test isempty(header)
        push!(header, SAM.MetaInfo("@HD\tVN:1.0\tSO:coordinate"))
        @test !isempty(header)
        @test length(header) == 1
        push!(header, SAM.MetaInfo("@CO\tsome comment"))
        @test length(header) == 2
        @test isa(collect(header), Vector{SAM.MetaInfo})
    end

    @testset "Record" begin
        record = SAM.Record()
        @test !isfilled(record)
        @test !SAM.ismapped(record)
        @test repr(record) == "XAM.SAM.Record: <not filled>"
        @test_throws ArgumentError SAM.flag(record)

        record = SAM.Record("r001\t99\tchr1\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*")
        @test isfilled(record)
        @test occursin(r"^XAM.SAM.Record:\n", repr(record))
        @test SAM.ismapped(record)
        @test SAM.isprimary(record)
        @test SAM.hastempname(record)
        @test SAM.tempname(record) == "r001"
        @test SAM.hasflag(record)
        @test SAM.flag(record) === UInt16(99)
        @test SAM.hasrefname(record)
        @test SAM.refname(record) == "chr1"
        @test SAM.hasposition(record)
        @test SAM.position(record) === 7
        @test SAM.hasmappingquality(record)
        @test SAM.mappingquality(record) === UInt8(30)
        @test SAM.hascigar(record)
        @test SAM.cigar(record) == "8M2I4M1D3M"
        @test SAM.hasnextrefname(record)
        @test SAM.nextrefname(record) == "="
        @test SAM.hasnextposition(record)
        @test SAM.nextposition(record) === 37
        @test SAM.hastemplength(record)
        @test SAM.templength(record) === 39
        @test SAM.hassequence(record)
        @test SAM.sequence(record) == dna"TTAGATAAAGGATACTG"
        @test !SAM.hasquality(record)
        @test_throws MissingFieldException SAM.quality(record)
    end

    @testset "Reader" begin
        reader = open(SAM.Reader, joinpath(samdir, "ce#1.sam"))
        @test isa(reader, SAM.Reader)
        @test eltype(reader) === SAM.Record

        # header
        h = header(reader)
        @test string.(findall(header(reader), "SQ")) == ["@SQ\tSN:CHROMOSOME_I\tLN:1009800"]

        # first record
        record = SAM.Record()
        read!(reader, record)
        @test SAM.ismapped(record)
        @test SAM.refname(record) == "CHROMOSOME_I"
        @test SAM.position(record) == leftposition(record) == 2
        @test SAM.rightposition(record) == rightposition(record) == 102
        @test SAM.tempname(record) == seqname(record) == "SRR065390.14978392"
        @test SAM.sequence(record) == sequence(record) == dna"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"
        @test SAM.sequence(String, record)          ==    "CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"
        @test SAM.seqlength(record) == 100
        @test SAM.quality(record)         == (b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" .- 33)
        @test SAM.quality(String, record) ==   "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        @test SAM.flag(record) == 16
        @test SAM.cigar(record) == "27M1D73M"
        @test SAM.alignment(record) == Alignment([
            AlignmentAnchor(  0,   1, OP_START),
            AlignmentAnchor( 27,  28, OP_MATCH),
            AlignmentAnchor( 27,  29, OP_DELETE),
            AlignmentAnchor(100, 102, OP_MATCH)])
        @test record["XG"] == 1
        @test record["XM"] == 5
        @test record["XN"] == 0
        @test record["XO"] == 1
        @test record["AS"] == -18
        @test record["XS"] == -18
        @test record["YT"] == "UU"
        @test eof(reader)
        close(reader)

        # rightposition (also implicitly alignlength)
        records = collect(open(SAM.Reader, joinpath(samdir, "ce#5b.sam")))
        @test SAM.rightposition(records[6]) == rightposition(records[6]) == 83

        # iterator
        @test length(collect(open(SAM.Reader, joinpath(samdir, "ce#1.sam")))) == 1
        @test length(collect(open(SAM.Reader, joinpath(samdir, "ce#2.sam")))) == 2

        # IOStream
        @test length(collect(SAM.Reader(open(joinpath(samdir, "ce#1.sam"))))) == 1
        @test length(collect(SAM.Reader(open(joinpath(samdir, "ce#2.sam"))))) == 2
    end

    @testset "Round trip" begin
        function compare_records(xs, ys)
            if length(xs) != length(ys)
                return false
            end
            for (x, y) in zip(xs, ys)
                if x.data[x.filled] != y.data[y.filled]
                    return false
                end
            end
            return true
        end
        for specimen in list_valid_specimens("SAM")
            filepath = joinpath(samdir, filename(specimen))
            mktemp() do path, io
                # copy
                reader = open(SAM.Reader, filepath)

                header_original = header(reader)

                writer = SAM.Writer(io, header_original)

                records = SAM.Record[]
                for record in reader
                    push!(records, record)
                    write(writer, record)
                end

                close(reader)
                close(writer)

                reader = open(SAM.Reader, path)

                @test header(reader) == header_original
                @test compare_records(collect(reader), records)

                close(reader)

            end
        end
    end
end

@testset "BAM" begin
    bamdir = path_of_format("BAM")

    @testset "AuxData" begin
        auxdata = BAM.AuxData(UInt8[])
        @test isempty(auxdata)

        buf = IOBuffer()
        write(buf, "NM", UInt8('s'), Int16(1))
        auxdata = BAM.AuxData(take!(buf))
        @test length(auxdata) == 1
        @test auxdata["NM"] === Int16(1)
        @test collect(auxdata) == ["NM" => Int16(1)]

        buf = IOBuffer()
        write(buf, "AS", UInt8('c'), Int8(-18))
        write(buf, "NM", UInt8('s'), Int16(1))
        write(buf, "XA", UInt8('f'), Float32(3.14))
        write(buf, "XB", UInt8('Z'), "some text\0")
        write(buf, "XC", UInt8('B'), UInt8('i'), Int32(3), Int32[10, -5, 8])
        auxdata = BAM.AuxData(take!(buf))
        @test length(auxdata) == 5
        @test auxdata["AS"] === Int8(-18)
        @test auxdata["NM"] === Int16(1)
        @test auxdata["XA"] === Float32(3.14)
        @test auxdata["XB"] == "some text"
        @test auxdata["XC"] == Int32[10, -5, 8]
        @test convert(Dict{String,Any}, auxdata) == Dict(
            "AS" => Int8(-18),
            "NM" => Int16(1),
            "XA" => Float32(3.14),
            "XB" => "some text",
            "XC" => Int32[10, -5, 8])
    end

    @testset "Record" begin
        record = BAM.Record()
        @test !isfilled(record)
        @test repr(record) == "XAM.BAM.Record: <not filled>"
        @test_throws ArgumentError BAM.flag(record)
    end

    @testset "Reader" begin
        reader = open(BAM.Reader, joinpath(bamdir, "ce#1.bam"))
        @test isa(reader, BAM.Reader)
        @test eltype(reader) === BAM.Record
        @test startswith(repr(reader), "XAM.BAM.Reader{IOStream}:")

        # header
        h = header(reader)
        @test isa(h, SAM.Header)

        # first record
        record = BAM.Record()
        read!(reader, record)
        @test BAM.ismapped(record)
        @test BAM.isprimary(record)
        @test ! BAM.ispositivestrand(record)
        @test BAM.refname(record) == "CHROMOSOME_I"
        @test BAM.refid(record) === 1
        @test BAM.hasnextrefid(record)
        @test BAM.nextrefid(record) === 0
        @test BAM.hasposition(record) === hasleftposition(record) === true
        @test BAM.position(record) === leftposition(record) === 2
        @test BAM.hasnextposition(record)
        @test BAM.nextposition(record) === 0
        @test rightposition(record) == 102
        @test BAM.hastempname(record) === hasseqname(record) === true
        @test BAM.tempname(record) == seqname(record) == "SRR065390.14978392"
        @test BAM.hassequence(record) === hassequence(record) === true
        @test BAM.sequence(record) == sequence(record) == dna"""
        CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCT
        AAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
        """
        @test BAM.seqlength(record) === 100
        @test BAM.hasquality(record)
        @test eltype(BAM.quality(record)) == UInt8
        @test BAM.quality(record) == [Int(x) - 33 for x in "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]
        @test BAM.flag(record) === UInt16(16)
        @test BAM.cigar(record) == "27M1D73M"
        @test BAM.alignment(record) == Alignment([
            AlignmentAnchor(  0,   1, OP_START),
            AlignmentAnchor( 27,  28, OP_MATCH),
            AlignmentAnchor( 27,  29, OP_DELETE),
            AlignmentAnchor(100, 102, OP_MATCH)])
        @test record["XG"] == 1
        @test record["XM"] == 5
        @test record["XN"] == 0
        @test record["XO"] == 1
        @test record["AS"] == -18
        @test record["XS"] == -18
        @test record["YT"] == "UU"
        @test keys(record) == ["XG","XM","XN","XO","AS","XS","YT"]
        @test values(record) == [1, 5, 0, 1, -18, -18, "UU"]
        @test eof(reader)
        close(reader)

        # Test conversion from byte array to record
        dsize = BAM.data_size(record)
        array = Vector{UInt8}(undef, BAM.FIXED_FIELDS_BYTES + dsize)
        GC.@preserve array record begin
            ptr = Ptr{UInt8}(pointer_from_objref(record))
            unsafe_copyto!(pointer(array), ptr, BAM.FIXED_FIELDS_BYTES)
            unsafe_copyto!(array, BAM.FIXED_FIELDS_BYTES + 1, record.data, 1, dsize)
        end
        new_record = convert(BAM.Record, array)
        @test record.bin_mq_nl == new_record.bin_mq_nl
        @test record.block_size == new_record.block_size
        @test record.flag_nc == new_record.flag_nc
        @test record.l_seq == new_record.l_seq
        @test record.next_refid == new_record.next_refid
        @test record.next_pos == new_record.next_pos
        @test record.refid == new_record.refid
        @test record.pos == new_record.pos
        @test record.tlen == new_record.tlen
        @test record.data == new_record.data

        # rightposition (also implicitly alignlength)
        records = collect(open(BAM.Reader, joinpath(bamdir, "ce#5b.bam")))
        @test BAM.rightposition(records[6]) == rightposition(records[6]) == 83

        # iterator
        @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#1.bam")))) == 1
        @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#2.bam")))) == 2

        # IOStream
        @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#1.bam"))))) == 1
        @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#2.bam"))))) == 2
    end

    @testset "Read long CIGARs" begin
        function get_cigar_lens(rec::BAM.Record)
            cigar_ops, cigar_n = BAM.cigar_rle(rec)
            field_ops, field_n = BAM.cigar_rle(rec, false)
            cigar_l = length(cigar_ops)
            field_l = length(field_ops)
            return cigar_l, field_l
        end

        function check_cigar_vs_field(rec::BAM.Record)
            cigar = BAM.cigar(rec)
            field = BAM.cigar(rec, false)
            cigar_l, field_l = get_cigar_lens(rec)
            return cigar != field && cigar_l != field_l
        end

        function check_cigar_lens(rec::BAM.Record, field_len, cigar_len)
            cigar_l, field_l = get_cigar_lens(rec)
            return cigar_l == cigar_len && field_l == field_len
        end

        reader = open(BAM.Reader, joinpath(bamdir, "cigar-64k.bam"))
        rec = BAM.Record()
        read!(reader, rec)
        @test !check_cigar_vs_field(rec)
        read!(reader, rec)
        @test check_cigar_vs_field(rec)
        @test check_cigar_lens(rec, 2, 72091)
    end

    function compare_records(xs, ys)
        if length(xs) != length(ys)
            return false
        end
        for (x, y) in zip(xs, ys)
            if !(
                x.block_size == y.block_size &&
                x.refid      == y.refid &&
                x.pos        == y.pos &&
                x.bin_mq_nl  == y.bin_mq_nl &&
                x.flag_nc    == y.flag_nc &&
                x.l_seq      == y.l_seq &&
                x.next_refid == y.next_refid &&
                x.next_pos   == y.next_pos &&
                x.tlen       == y.tlen &&
                x.data[1:BAM.data_size(x)] == y.data[1:BAM.data_size(y)])
                return false
            end
        end
        return true
    end

    @testset "Round trip" begin
        for specimen in list_valid_specimens("BAM")
            filepath = joinpath(bamdir, filename(specimen))
            mktemp() do path, _
                # copy
                if hastags(specimen) && in("bai", tags(specimen))
                    reader = open(BAM.Reader, filepath, index=filepath * ".bai")
                else
                    reader = open(BAM.Reader, filepath)
                end

                header_original = header(reader)

                writer = BAM.Writer(BGZFStream(path, "w"), BAM.header(reader, fillSQ=isempty(findall(header(reader), "SQ"))))

                records = BAM.Record[]
                for record in reader
                    push!(records, record)
                    write(writer, record)
                end
                close(reader)
                close(writer)

                reader = open(BAM.Reader, path)

                @test header(reader) == header_original
                @test compare_records(collect(reader), records)

                close(reader)

            end
        end
    end

    @testset "Random access" begin
        filepath = joinpath(bamdir, "GSE25840_GSM424320_GM06985_gencode_spliced.head.bam")
        reader = open(BAM.Reader, filepath, index=filepath * ".bai")

        @test isa(eachoverlap(reader, "chr1", 1:100), BAM.OverlapIterator)
        @test isa(eachoverlap(reader, GenomicFeatures.Interval("chr1", 1, 100)), BAM.OverlapIterator)

        # expected values are counted using samtools
        for (refname, interval, expected) in [
                ("chr1", 1_000:10000,      21),
                ("chr1", 8_000:10000,      20),
                ("chr1", 766_000:800_000, 142),
                ("chr1", 786_000:800_000, 1),
                ("chr1", 796_000:800_000, 0)]
            intsect = eachoverlap(reader, refname, interval)
            @test eltype(intsect) == BAM.Record
            @test count(_ -> true, intsect) == expected
            # check that the intersection iterator is stateless
            @test count(_ -> true, intsect) == expected
        end

        # randomized tests
        for n in 1:50
            refindex = 1
            refname = "chr1"
            range = randrange(1:1_000_000)
            seekstart(reader)
            # linear scan
            expected = filter(collect(reader)) do record
                BAM.compare_intervals(record, (refindex, range)) == 0
            end
            # indexed scan
            actual = collect(eachoverlap(reader, refname, range))
            @test compare_records(actual, expected)
        end
        close(reader)

        filepath = joinpath(bamdir, "R_12h_D06.uniq.q40.bam")
        reader = open(BAM.Reader, filepath, index=filepath * ".bai")
        @test isempty(collect(eachoverlap(reader, "chr19", 5823708:5846478)))
        close(reader)
    end
end
