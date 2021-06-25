function quality_string(quals::Vector{UInt8})
    characters = Vector{Char}(undef, length(quals))
    for i in eachindex(quals)
        @inbounds qual = quals[i]
        if qual < 10
            char = ' '
        elseif qual < 15
            char = '▁'
        elseif qual < 20
            char = '▂'
        elseif qual < 25
            char = '▃'
        elseif qual < 30
            char = '▄'
        elseif qual < 35
            char = '▆'
        elseif qual < 40
            char = '▇'
        elseif qual < 255
            char = '█'
        else
            char = '?'
        end
        @inbounds characters[i] = char
    end
    return join(characters)
end

function compact_string(sequence, width)
    if length(sequence) <= width
        return string(sequence)
    else
        half = div(width - 1, 2)
        return string(sequence[1:half]) * '…' * string(sequence[end-half:end])
    end
end
