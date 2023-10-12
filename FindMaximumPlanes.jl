########################################################################
########################################################################

function _FindMaximumPlane(Pairs::Array{Tuple,1},label::Array{Dict{Any,Any},1},df,last::Int64)
    
    lab1 = df[1:last,2]
    lab2 = df[last+1:end,2]
    
    
    for p in Pairs
        if label[1]
        end
    end
end

########################################################################
########################################################################

function FindMaximumPlane(Pairs::Array{Any,1},label::Array{Dict{Any,Any},1},CSVpathName)
    
    nPlane = length(label)
    maxZ = div(nPlane,2)+1
    
    df, zMin, zMax, first, last = FindPlanes(CSVpathName)
    
    for z in 1:div(nPlane,2)
        if !isempty(label[maxZ-z]) && !isempty(label[maxZ-z+1])
            _FindMaximumIntensityPlane(Pairs[maxZ-z],label[maxZ-z:maxZ-z+1],df[first[maxZ-z]:last[maxZ-z+1],:],last[maxZ-z])
        end
        if !isempty(label[maxZ+z-1]) && !isempty(label[maxZ+z])
            _FindMaximumIntensityPlane(Pairs[maxZ+z-1],label[maxZ+z-1:maxZ+z],df[first[maxZ+z-1]:last[maxZ+z],:],last[maxZ+z-1])
        end
    end
    
    return
    
end

########################################################################
########################################################################