########################################################################
########################################################################

using TiffImages, Images

########################################################################
########################################################################

function Shift!(label::Dict{Any,Any},currentKey::Int64,nextKey::Int64)
    
    shiftAmplitude = nextKey-currentKey
    
    for k in nextKey:maximum(keys(label))
        
        if haskey(label,k)
            label[k-shiftAmplitude] = label[k]
            delete!(label,k)
        end
    end
    
    return
end

########################################################################
########################################################################

function Relabelling!(label::Dict{Any,Any}) # It doesn't work properly
    
    orderedKeys = sort([k for k in keys(label)])
    nextKey = 1
    currentKey = 1
    currentMax = maximum(keys(label))
        
    while currentKey<=currentMax
        if !haskey(label,currentKey)
            
            Shift!(label,currentKey,orderedKeys[nextKey])
            orderedKeys = sort([k for k in keys(label)])
            currentMax = maximum(keys(label))
                
        end
        currentKey += 1
        nextKey += 1
    end
    
    return
end

########################################################################
########################################################################

function NewDroplet!(coord::Array{Int64,1},yTemp::Int64,label::Dict{Any,Any},label2::Dict{Any,Any})
    
    x,y,z = coord
    
    if !isempty(keys(label))
        
        # Even if we didn't find a droplet, all the processed pixels in the line can be assigned to the same
        l = maximum(values(label))+1
        label2[l] = []
        for Y in y:yTemp
            label[[x,Y,z]] = l
            push!(label2[l],[x,Y,z])
        end
    else
        label2[1] = []
        for Y in y:yTemp
            label[[x,Y,z]] = 1
            push!(label2[1],[x,Y,z])
        end
    end
    
    return
    
end

########################################################################
########################################################################

function Queue(coord::Array{Int64,1},dense,label::Dict{Any,Any},label2::Dict{Any,Any})
    
    x,y,z = coord
    xy = size(dense)
    flag = -1
    yTemp = y
    
    while flag < 0
        yTemp+=1
        if dense[x,yTemp] == 1
            # in this case we found the label, not only for the starting pixel but for all the pixels in the line with Y in between the first and last analysed
            if dense[x-1,yTemp] == 1
                for Y in y:yTemp
                    label[[x,Y,z]] = label[[x-1,yTemp,z]]
                    push!(label2[label[[x,Y,z]]],[x,Y,z])
                end
                flag = 0
            end
        else
            # we finished the line in the droplet without a label found -> new droplet
            flag=1
        end
        # we finished the line in the image without a label found -> new droplet
        if yTemp == xy[2]
            flag = 1
        end
    end
    
    return flag, yTemp
    
end

########################################################################
########################################################################

function MergeDroplets!(coord::Array{Int64,1},label::Dict{Any,Any},label2::Dict{Any,Any})
    
    x,y,z = coord
    l = min(label[[x,y-1,z]], label[[x-1,y,z]])
    L = max(label[[x,y-1,z]], label[[x-1,y,z]])
    
    # merge the droplets
    for v in label2[L]
        push!(label2[l],v)
        delete!(label,v)
        label[v] = l
    end
    push!(label2[l],[x,y,z])
    label[[x,y,z]] = l
    delete!(label2,L)
    
    return
    
end

########################################################################
########################################################################

function labellingLine!(x::Int64,dense,z::Int64,label::Dict{Any,Any},label2::Dict{Any,Any})
    
    xy = size(dense)
    
    for y in 1:xy[2]
        if dense[x,y] == 1

            # divide the droplets in the plane
            # look for bright pixels in the left upper region
            if x>1 && dense[x-1,y] == 1
                if y>1 && dense[x,y-1] == 1 && label[[x,y-1,z]] != label[[x-1,y,z]]
                    MergeDroplets!([x,y,z],label,label2)
                else
                    label[[x,y,z]]=label[[x-1,y,z]]
                    push!(label2[label[[x,y,z]]],[x,y,z])
                end
            elseif y>1 && dense[x,y-1] == 1
                if x>1 && dense[x-1,y] == 1 && label[[x-1,y,z]] != label[[x,y-1,z]]
                    MergeDroplets!([x,y,z],label,label2)
                else
                    label[[x,y,z]]=label[[x,y-1,z]]
                    push!(label2[label[[x,y,z]]],[x,y,z])
                end

            # It could be that a part of this droplets is already present in the upper right region, if this is the case the new pixel doesn't belong to a new droplets but we need a queue to find the right one. In this case we'll process all the bright pixel on the right in the same line and we'll look at the corresponding pixel above, if it's found we have the label otherwise, if we find a dark pixel before a label we are in a new droplet.
            elseif x>1 && y<xy[2] && (dense[x,y+1] == 1)
                flag, yTemp = Queue([x,y,z],dense,label,label2)
                # if we exit with 0 or 1 all the pixels in the line with y<yTemp belong to the same droplet, so let skip them
                if flag == 0
                    y = yTemp
                elseif flag == 1
                    NewDroplet!([x,y,z],yTemp,label,label2)
                    y = yTemp
                end
            else
                NewDroplet!([x,y,z],y,label,label2)
            end
        end
    end
    
    return
    
end
    
########################################################################
########################################################################

function _Labelling!(dense,z::Int64,label::Dict{Any,Any},label2::Dict{Any,Any})
    
    xy = size(dense)
    
    for x in 1:xy[1]
        labellingLine!(x,dense,z,label,label2)
    end
    
    return
    
end

########################################################################
########################################################################

function Labelling(dense,intensity,nPlane::Int64)
    
    label = Array{Dict{Any,Any},1}(undef,nPlane+1)
    label2 = Array{Dict{Any,Any},1}(undef,nPlane+1)
    
    # label droplets in maxZ ± nPlane/2
    maxIntZ = findmax(intensity)[2]
    minZ = maxIntZ-div(nPlane,2)
    maxZ = maxIntZ+div(nPlane,2)
    if minZ < 1
        minZ = 1
    end
    if maxZ > length(intensity)
        maxZ = length(intensity)
    end
    println("Maximum intensity plane ",maxIntZ)
    
    i = 0
    for z in minZ:maxZ
        i += 1
        label[i] = Dict()
        label2[i] = Dict()
        if intensity[z]>0
            println("Labelling plane ",z)
            _Labelling!(dense[:,:,z],z,label[i],label2[i])

            # Check if some index is unused and relabel droplets in the plane
            if maximum(keys(label2[i]))>length(keys(label2[i]))

                Relabelling!(label2[i])
            end
        else
            println("Zero intensity plane ",z)
        end
    end
    
    return label2
    
end

########################################################################
########################################################################