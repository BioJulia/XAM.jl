# SAM File Format
# ===============

module SAM

using BioGenerics

import Automa
import Automa.RegExp: @re_str
import BioAlignments
import BioGenerics.Exceptions: missingerror
import BioGenerics.RecordHelper: unsafe_parse_decimal
import BioGenerics: isfilled, header
import BioSequences
import BufferedStreams
using Printf: @sprintf

include("flags.jl")
include("metainfo.jl")
include("record.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
