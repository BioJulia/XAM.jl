@testset "Issues" begin

    # https://github.com/BioJulia/XAM.jl/issues/31
    path_bam = joinpath(path_of_format("BAM"), "SRR7993829_1.100K.forward.bam")

    open(BAM.Reader, path_bam, index = path_bam * ".bai") do reader

        @test count(overlap -> true, eachoverlap(reader, "JH584304.1", 51000:51200)) == 0
        @test count(overlap -> true, eachoverlap(reader, "JH584304.1", 51000:51715)) == 1

    end

end
