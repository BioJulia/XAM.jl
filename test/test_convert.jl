using Test, XAM

@testset "convert SAM to BAM" begin

    samdir = path_of_format("SAM")

    function test_auxdata(record, bam_record)
        aux_sam = SAM.auxdata(record)
        aux_bam = BAM.auxdata(bam_record)

        @test Set(keys(aux_sam)) == Set(keys(aux_bam))
        @test all(aux_sam[k] == aux_bam[k] for k in keys(aux_sam))
    end

    function test_record(record, bam_record)
        @test SAM.position(record) == BAM.position(bam_record)
        @test SAM.flags(record) == BAM.flags(bam_record)
        @test SAM.cigar(record) == BAM.cigar(bam_record)
        @test SAM.sequence(record) == BAM.sequence(bam_record)
        if SAM.hasquality(record)
            @test SAM.quality(record) == BAM.quality(bam_record)
        end
        test_auxdata(record, bam_record)
    end

    # internal methods
    @test XAM.Conversion.reg2bin(-1,0) == 4680

    # good alignment with missing qualities
    sam = "seq1\t81\tPhiX\t1051\t60\t70M\t=\t1821\t702\tTCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC\t*\tNM:i:0\tMD:Z:70\tMC:Z:70M\tAS:i:70\tXS:i:0"
    record = XAM.SAM.Record(sam)

    bam_record = BAM.Record(record)
    test_record(record, bam_record)
    @test BAM.quality(bam_record) == fill(0xff, BAM.seqlength(bam_record))

    # good alignment with qualities
    sam = "seq1\t81\tPhiX\t1051\t60\t70M\t=\t1821\t702\tTCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC\tADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tNM:i:0\tMD:Z:70\tMC:Z:70M\tAS:i:70\tXS:i:0"
    record = XAM.SAM.Record(sam)
    bam_record = BAM.Record(record)
    test_record(record, bam_record)    
    @test BAM.quality(bam_record) == BAM.quality(bam_record)

    # with an insertion
    sam = "seq1\t16\tPhiX\t1051\t60\t38M6I32M\t*\t0\t0\tTCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCCCCCCCATTTCAACTACTCCGGTTATCGCTGGCGAC\t*\tNM:i:6\tMD:Z:70\tAS:i:58\tXS:i:0"
    
    record = XAM.SAM.Record(sam)
    bam_record = BAM.Record(record)
    test_record(record, bam_record)

    #################
    ## test all SAM files in FormatSpecimens
    function test_sam_file(samdir, file)

        records = collect(open(SAM.Reader, joinpath(samdir, file)))
        
        for record in records
            # record with more than 65535 operations
            if SAM.tempname(record) == "03d98240-9a1d-4a9a-ae95-2fa8e2f1c945_Basecall_1D_template"
                @test_throws ErrorException BAM.Record(record)
                continue
            end
            
            bam_record = BAM.Record(record)
            test_record(record, bam_record)
        end
    end

    for file in ["ce#1.sam","ce#2.sam","ce#5.sam","ce#5b.sam","ce#unmap.sam","ce#unmap1.sam","ce#unmap2.sam","cigar-64k.sam","ex1.sam","ex1_header.sam","sam1.sam","sam2.sam","xx#blank.sam","xx#minimal.sam"]
        test_sam_file(samdir, file)
    end

    #records = collect(open(SAM.Reader, joinpath(samdir, "ex1.sam")))

    #################
    # test with a BAM

    function test_bam(bam_file, sam_file)
        bamdir = path_of_format("BAM")
        reader = open(BAM.Reader, joinpath(bamdir, bam_file))
        h = BAM.header(reader)
        bam_records = [record for record in reader] 

        reader = open(SAM.Reader, joinpath(samdir, sam_file))
        sam_records = [record for record in reader] 

        for (sam_record, bam_record) in zip(sam_records, bam_records)
            if SAM.tempname(sam_record) == "03d98240-9a1d-4a9a-ae95-2fa8e2f1c945_Basecall_1D_template"
                continue
            end
            #@info SAM.tempname(sam_record)
            bam_record2 = convert(BAM.Record, sam_record; header=h)
            @test bam_record == bam_record2
        end

    end

    test_bam("ce#1.bam", "ce#1.sam")
    test_bam("ce#2.bam", "ce#2.sam")
    test_bam("ce#5.bam", "ce#5.sam")
    test_bam("ce#unmap.bam", "ce#unmap.sam")
    test_bam("ce#unmap1.bam", "ce#unmap1.sam")
    test_bam("ce#unmap2.bam", "ce#unmap2.sam")
    #test_bam("cigar-64k.bam", "cigar-64k.sam") # this one as some issues
    
end