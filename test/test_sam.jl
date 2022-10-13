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
            AlignmentAnchor(  0,   1,   0, OP_START),
            AlignmentAnchor( 27,  28,  27, OP_MATCH),
            AlignmentAnchor( 27,  29,  28, OP_DELETE),
            AlignmentAnchor(100, 102, 101, OP_MATCH)])
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

    @testset "In-Place-Reading Pattern" begin

        file_sam = joinpath(samdir, "ce#5b.sam")

        records = open(collect, SAM.Reader, file_sam)

        reader = open(SAM.Reader, file_sam)
        record = SAM.Record()
        i = 0
        while !eof(reader)
            empty!(record) # Reset the record.
            read!(reader, record)

            i = i + 1

            @test records[i] == record

        end

        close(reader)

        # Test blank file.
        file_sam = joinpath(samdir, "xx#blank.sam")

        records = open(collect, SAM.Reader, file_sam)
        @test records == []

    end
end
