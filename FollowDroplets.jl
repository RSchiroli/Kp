########################################################################
########################################################################

function Distance(point1,point2)
    return sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2)
end

########################################################################
########################################################################

# The idea behind this function is that if two droplets in two planes are the same, a rand pixel in the first and another in the second should be at a distance comparable with the diameter of the droplet itself
function NotTooFar(droplet::Array{Any,1},labelPlane2::Dict{Any,Any})
    
    sizeDroplet = length(droplet)
    dropletsToCheck = []
    p1 = rand(droplet)
    
    for k2 in keys(labelPlane2)
        p2 = rand(labelPlane2[k2])
        if Distance(p1,p2) < sqrt(max(sizeDroplet,length(labelPlane2[k2])))*100 # <--- TO BE IMPROVED
            push!(dropletsToCheck,k2)
        end
    end
    
    return dropletsToCheck
    
end

########################################################################
########################################################################

function _FollowDroplets(labelPlane1::Dict{Any,Any},labelPlane2::Dict{Any,Any})
    
    pairs = fill((-1,-1),length(keys(labelPlane1)))
    
    for k in keys(labelPlane1)
        found = 0
        dropletsToCheck = NotTooFar(labelPlane1[k],labelPlane2)
        for pixel in labelPlane1[k]
            found == 0 || break
            newCoord = [pixel[1],pixel[2],pixel[3]+1]
            for k2 in dropletsToCheck
                if newCoord in labelPlane2[k2]
                    pairs[k] = (k,k2)
                    found = 1
                end
            end
        end
        if found == 0
            println("droplet ",k," not found")
        end
    end
    
    
    
    return pairs# vector of tuple with pair of indices (first_label, last_label) ?
end

########################################################################
########################################################################

function FollowDroplets(label::Array{Dict{Any,Any},1})
    
    nPlanes = length(label)
    nPlanes>1 || @error "multiple planes needed to follow droplets"
    
    Pairs = []
    
    for i in 1:nPlanes-1
        push!(Pairs,_FollowDroplets(label[i],label[i+1]))
            
    end
    
    return Pairs
    
end

########################################################################
########################################################################