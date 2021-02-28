@testset "Cross Check Properties" begin

	Broadcast.broadcastable(x::XAM.BAM.Record) = Ref(x) #TODO: consider moving to XAM.jl.
	Broadcast.broadcastable(x::XAM.SAM.Record) = Ref(x) #TODO: consider moving to XAM.jl

	function crosscheck(bam::BAM.Record, sam::SAM.Record, property::Symbol)

		bam_property = getproperty(XAM.BAM, property)
		sam_property = getproperty(XAM.SAM, property)

		if bam_property(bam) != sam_property(sam)
			@warn "$property result is not the same" bam_property(bam) sam_property(sam)
			return false
		end

		return true
	end

	samdir = path_of_format("SAM")
	bamdir = path_of_format("BAM")

	filenames = [
		"ce#1",
		"ce#2",
		"ce#5",
		"ce#5b",
		"ce#unmap",
		"ce#unmap1",
		"ce#unmap2",
	]

	properties = [
		:position,# POS
		:tempname,# QNAME
		:mappingquality,# MAPQ
		:cigar, # CIGAR
		:flag, # FLAG
		:sequence, # SEQ
		:nextposition, # PNEXT
		:templength, # TLEN
	]

	for filename in filenames

		records_bam = collect(open(BAM.Reader, joinpath(bamdir, filename * ".bam")))
		records_sam = collect(open(SAM.Reader, joinpath(samdir, filename * ".sam")))

		for (bam, sam) in zip(records_bam, records_sam)
			@test all(crosscheck.(bam, sam, properties)) == true
		end

	end

end # testset Crosscheck
