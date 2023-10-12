########################################################################
########################################################################

include("ComputeKp.jl")
include("FollowDroplets.jl")
include("BlockStatistics.jl")
include("FindMaximumPlanes.jl")

########################################################################
########################################################################

function RunKpAnalysis(imgName,maskName;Gain::Int64=100,threshold::Int64=0,nPlane::Int64=0,returnAll::Bool=false)
    
    avgKp,stdKp,allKp,Sizes,label = CalcKp(imgName,maskName,Gain=Gain,threshold=threshold,nPlane=nPlane,returnAll=true)
    
    # Pairs = FollowDroplets(label)
    
    path = split(imgName,"/")
    name = path[end]
    fileName,_ = split(name,".")
    csvPath = ""
    for folder in path[1:end-1]
        csvPath *= folder*"/"
    end
    csvPath*=fileName*".csv"
    
    # result = BlockStats(csvPath)
    
    if returnAll==true
        return avgKp,stdKp,allKp,Sizes,label#Pairs,result
    end
    
    # return result
    
    return
    
end

########################################################################
########################################################################