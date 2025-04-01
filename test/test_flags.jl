@testset "Flags" begin

    flag_tests = [
        (XAM.ispaired, XAM.FLAG_PAIRED),
        (XAM.isproperpair, XAM.FLAG_PROPER_PAIR),
        (XAM.isunmapped, XAM.FLAG_UNMAPPED), 
        (XAM.ismapped, ~XAM.FLAG_UNMAPPED), 
        (XAM.isnextunmapped, XAM.FLAG_NEXT_UNMAPPED), 
        (XAM.isnextmapped, ~XAM.FLAG_NEXT_UNMAPPED), 
        (XAM.isreversecomplemented, XAM.FLAG_REVERSE), 
        (XAM.isforwardstrand, ~XAM.FLAG_REVERSE), 
        (XAM.ispositivestrand, ~XAM.FLAG_REVERSE), 
        (XAM.isreversestrand, XAM.FLAG_REVERSE), 
        (XAM.isnegativestrand, XAM.FLAG_REVERSE), 
        (XAM.isnextreversecomplemented, XAM.FLAG_NEXT_REVERSE), 
        (XAM.isfirstsegment, XAM.FLAG_FIRST_SEGMENT), 
        (XAM.isread1, XAM.FLAG_FIRST_SEGMENT), 
        (XAM.islastsegment, XAM.FLAG_LAST_SEGMENT), 
        (XAM.isread2, XAM.FLAG_LAST_SEGMENT), 
        (XAM.issecondaryalignment, XAM.FLAG_SECONDARY), 
        (XAM.isqcfail, XAM.FLAG_QCFAIL), 
        (XAM.isduplicate, XAM.FLAG_DUPLICATE), 
        (XAM.issupplementaryalignment, XAM.FLAG_SUPPLEMENTARY), 
        (XAM.isprimaryalignment, ~0x900), 
    ]

    foreach(flag_tests) do (func, flag)
        # SAM Record
        io = IOBuffer()
        write(io, string(Int(flag)))
        bytes = take!(io)

        sam_record = XAM.SAM.Record()
        sam_record.data = bytes
        sam_record.flags = sam_record.filled = 1:length(bytes)

        @test flag == XAM.flags(sam_record)
        @test func(sam_record)

        # BAM record
        bam_record = XAM.BAM.Record(sam_record)

        @test typeof(bam_record) == XAM.BAM.Record
        @test flag == XAM.flags(bam_record)
        @test func(bam_record)
    end

end