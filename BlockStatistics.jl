########################################################################
########################################################################

using CSV, Tables, DataFrames, PyPlot

########################################################################
########################################################################

mutable struct result_to_plot
    size::Array{Int64,1}
    Kp::Array{Float64,1}
    intensity::Array{Float64,1}
    intensityDiluteRel::Float64
    nDilute::Int64
end

########################################################################
########################################################################

function ResultToPlot()
    size = Array{Int64,1}(undef,0)
    Kp = Array{Float64,1}(undef,0)
    intensity = Array{Float64,1}(undef,0)
    intensityDiluteRel = 0.0
    nDilute = 0
    return result_to_plot(size,Kp,intensity,intensityDiluteRel,nDilute)
end

########################################################################
########################################################################

function _BlockStats(kp::Array{Float64,1},sizes::Array{Float64,1})
    
    avg = sum(sizes)/length(sizes)
    std = sqrt((sum((sizes.-avg).^2)/length(sizes)))
    
    
    AvgKpBlocks = [0.0,0.0,0.0,0.0,0.0]
    StdKpBlocks = [0.0,0.0,0.0,0.0,0.0]
    AvgSizeBlocks = [0.0,0.0,0.0,0.0,0.0]
    StdSizeBlocks = [0.0,0.0,0.0,0.0,0.0]
    tempKpBlocks = [[],[],[],[],[]]
    tempSizeBlocks = [[],[],[],[],[]]
    DropletsPerBlock = [0,0,0,0,0]
    
    
    for i in 1:length(kp)
        if sizes[i]<avg-3*std/2
            AvgKpBlocks[1] += kp[i]
            AvgSizeBlocks[1] += sizes[i]
            push!(tempKpBlocks[1],kp[i])
            push!(tempSizeBlocks[1],sizes[i])
        elseif avg-3*std/2<=sizes[i]<avg-std/2
            AvgKpBlocks[2] += kp[i]
            AvgSizeBlocks[2] += sizes[i]
            push!(tempKpBlocks[2],kp[i])
            push!(tempSizeBlocks[2],sizes[i])
        elseif avg-std/2<=sizes[i]<avg+std/2
            AvgKpBlocks[3] += kp[i]
            AvgSizeBlocks[3] += sizes[i]
            push!(tempKpBlocks[3],kp[i])
            push!(tempSizeBlocks[3],sizes[i])
        elseif avg+std/2<=sizes[i]<avg+3*std/2
            AvgKpBlocks[4] += kp[i]
            AvgSizeBlocks[4] += sizes[i]
            push!(tempKpBlocks[4],kp[i])
            push!(tempSizeBlocks[4],sizes[i])
        else
            AvgKpBlocks[5] += kp[i]
            AvgSizeBlocks[5] += sizes[i]
            push!(tempKpBlocks[5],kp[i])
            push!(tempSizeBlocks[5],sizes[i])
        end
    end
    
    for i in 1:5
        if !isempty(tempKpBlocks[i])
            sizeBlock = length(tempKpBlocks[i])
            AvgKpBlocks[i] /= sizeBlock
            StdKpBlocks[i] = sqrt((sum((tempKpBlocks[i].-AvgKpBlocks[i]).^2)/sizeBlock))
            AvgSizeBlocks[i] /= sizeBlock
            StdSizeBlocks[i] = sqrt((sum((tempSizeBlocks[i].-AvgSizeBlocks[i]).^2)/sizeBlock))
            DropletsPerBlock[i] = length(tempKpBlocks[i])
        else
            DropletsPerBlock[i] = 0
        end
        
    end
        
    return AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks, DropletsPerBlock
    
end

########################################################################
########################################################################

function FindPlanes(CSVpathName)
    
    df = CSV.read(CSVpathName,DataFrame)
    zMin = Int64(minimum(df[:,1]))
    zMax = Int64(maximum(df[:,1]))
    # first and last are indices of first and last element in the plane corresponding to their index
    first = [1]
    last = []
    z = zMin
    for i in 1:size(df,1) 
        if df[i,1] == z+1
            push!(last,i-1)
            push!(first,i)
            z=df[i,1]
        end
        if z == zMax
            push!(last,size(df,1))
            i = size(df,1)
        end
    end
    
    return df, zMin, zMax, first, last
end

########################################################################
########################################################################

function BlockStats(CSVpathName)
    
    df, zMin, zMax, first, last = FindPlanes(CSVpathName)
    
    results = []
    stats = []
    
    for z in 1:zMax-zMin+1
        AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks, DropletsPerBlock = _BlockStats(df[first[z]:last[z],4],df[first[z]:last[z],3])
        push!(results,[AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks, DropletsPerBlock])
        push!(stats,[df[first[z],6],last[z]-first[z]+1,df[first[z],7]])
    end
    
    # results = avgKpBlock, StdKpBlock, AvgSizeBlock, StdSizeBlock
    # stats = intensityDilutePhaseInThePlane, NumberOfDropletsInThePlane, NumberOfPixels in the dilute phase
    
    return results,stats
    
end

########################################################################
########################################################################

function FindMaximumIntensityPlanes(CSVpathName)
    
    df = CSV.read(CSVpathName,DataFrame)
    zMin = Int64(minimum(df[:,1]))
    zMax = Int64(maximum(df[:,1]))
            
    # If there is only one plane labelled everything is simpler
    if zMin==zMax
        return df, zMin, 1, size(df,1)
    else
        # first and last are indices of first and last element in the plane corresponding to their index
        first = [1]
        last = []
        intensity = []
        z = zMin
        for i in 1:size(df,1) 
            if df[i,1] == z+1
                push!(last,i-1)
                push!(first,i)
                z=df[i,1]
                push!(intensity,sum(df[first[end]:last[end],5]))
                #
            end
            if z == zMax
                push!(last,size(df,1))
                i = size(df,1)
            end
        end
        # maxIntZ is the maximum intensity plane's index between the labelled, the returned value is the plane's index in the image
        _,maxIntZ = findmax(intensity)
    end
    
    return df, maxIntZ+zMin-1, first[maxIntZ], last[maxIntZ]
end

########################################################################
########################################################################

function NamePath(CSVpathName)
    path = split(CSVpathName,"/")
    Name = ""
    for folder in path[1:end-1]
        Name = Name*folder*"/"
    end
    return Name
end

########################################################################
########################################################################

function MultipleImagesBlockStats(CSVpathNames,Title)
    Res = []
    Result = []
    AllSizes = Array{Float64,1}(undef,0)
    AllKp = Array{Float64,1}(undef,0)
    indexMaxInt = []
    
    for name in CSVpathNames
        res = ResultToPlot()
        result = BlockStats(name)
        df, zMaxInt, first, last = FindMaximumIntensityPlanes(name)
        push!(indexMaxInt,div(length(first),2)+1)
        res.size = Int64.(df[first:last,3])
        res.Kp = df[first:last,4]
        res.intensity = df[first:last,5]
        res.intensityDiluteRel = df[first,6]
        res.nDilute = Int64(df[first,7])
        push!(Res,res)
        push!(Result,result)
        AllSizes = vcat(AllSizes,res.size)
        AllKp = vcat(AllKp,res.Kp)
    end
    
    AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks, DropletsPerBlock = _BlockStats(AllKp,AllSizes)
    
    path = NamePath(CSVpathNames[1])
    mkdir(path*Title)
    
    plotName="_avgKp_images.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    Legend = []
    # x = avgSize, y = avgKp, yerr = stdKp
    for i in 1:length(Result)
        errorbar(Result[i][1][indexMaxInt[i]][3],Result[i][1][indexMaxInt[i]][1],Result[i][1][indexMaxInt[i]][2],fmt=".")
        push!(Legend,string(i))
        #
    end
    xlabel("avgSize")
    ylabel("avgKp")
    title(Title*" "*plotName)
    legend(Legend)
    savefig(NameToSavePlot)
    close()
        
    plotName="_Kp_vs_DiluteInt.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    Legend = []
    # x = avgKp, y = diluteInt
    for i in 1:length(Result)
        scatter(Result[i][1][indexMaxInt[i]][1][2:end],[Result[i][2][indexMaxInt[i]][1] for j in 1:4])
    end
    xlabel("avgKp")
    ylabel("Dilute Intensity")
    title(Title*" "*plotName)
    savefig(NameToSavePlot)
    close()
    
    plotName="_DiluteInt_vs_nDilute.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    Legend = []
    # x = diluteInt, y = nDilute
    for i in 1:length(Result)
        scatter(Result[i][2][indexMaxInt[i]][3],Result[i][2][indexMaxInt[i]][1])
    end
    xlabel("n Dilute")
    ylabel("Dilute Intensity")
    title(Title*" "*plotName)
    savefig(NameToSavePlot)
    close()
    
    plotName="_avgKp_All.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    errorbar(AvgSizeBlocks,AvgKpBlocks,StdKpBlocks,fmt=".")
    xlabel("avgSize")
    ylabel("avgKp")
    title(Title*" "*plotName)
    savefig(NameToSavePlot)
    close()
    
    plotName="_avgSize.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    errorbar(LinRange(1,length(AvgSizeBlocks),length(AvgSizeBlocks)),AvgSizeBlocks,StdSizeBlocks,fmt=".")
    xlabel("n Block")
    ylabel("avgSize")
    title(Title*" "*plotName)
    savefig(NameToSavePlot)
    close()
    
    plotName="_BlockSize.svg"
    NameToSavePlot = path*Title*"/"*Title*plotName
    scatter(AvgSizeBlocks,DropletsPerBlock)
    xlabel("avgSize")
    ylabel("# droplets")
    title(Title*" "*plotName)
    savefig(NameToSavePlot)
    close()

    return AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks
    
end

########################################################################
########################################################################

function _Stats(kp::Array{Float64,1},sizes::Array{Float64,1})
    
    avgKp = sum(kp)/length(kp)
    stdKp = sqrt((sum((kp.-avgKp).^2)/length(kp)))
    
    avgSize = sum(sizes)/length(sizes)
    stdSize = sqrt((sum((sizes.-avgSize).^2)/length(sizes)))
            
    return avgKp, stdKp, avgSize, stdSize
    
end

########################################################################
########################################################################

function Stats(CSVpathName)
    
    df, zMin, zMax, first, last = FindPlanes(CSVpathName)
    
    results = []
    stats = []
    
    for z in 1:zMax-zMin+1
        AvgKp, StdKp, AvgSize, StdSize = _Stats(df[first[z]:last[z],4],df[first[z]:last[z],3])
        push!(results,[AvgKp, StdKp, AvgSize, StdSize])
#        push!(stats,[df[first[z],6],last[z]-first[z]+1,df[first[z],7]])
    end
       
    return results#,stats
    
end

########################################################################
########################################################################

function MultipleImagesStats(CSVpathNames,Title)
    
    Res = fill(0.0,length(CSVpathNames)+3,4)
    #Result = []
    AllSizes = Array{Float64,1}(undef,0)
    AllKp = Array{Float64,1}(undef,0)
    indexMaxInt = []
    
    i=0
    for name in CSVpathNames
        
        i+=1
        result = Stats(name)
     #   push!(Result,result)
        
        df, zMaxInt, first, last = FindMaximumIntensityPlanes(name)
        push!(indexMaxInt,div(length(first),2)+1)
        size = Int64.(df[first:last,3])
        Kp = df[first:last,4]
        
        Res[i,1] = result[indexMaxInt[end]][1]
        Res[i,2] = result[indexMaxInt[end]][2]
        Res[i,3] = result[indexMaxInt[end]][3]
        Res[i,4] = result[indexMaxInt[end]][4]
        
        AllSizes = vcat(AllSizes,size)
        AllKp = vcat(AllKp,Kp)
        
    end
    
    AvgKp, StdKp, AvgSize, StdSize = _Stats(AllKp,AllSizes)
    
    Res[end-1,1] = AvgKp
    Res[end-1,2] = StdKp
    Res[end-1,3] = AvgSize
    Res[end-1,4] = StdSize
    
    Res[end,1] = sum([Res[i,1] for i in 1:size(Res,1)])/size(Res,1)
    Res[end,2] = sqrt((sum([Res[i,1] for i in 1:size(Res,1)].-Res[end,1]).^2)/size(Res,1))
    Res[end,3] = sum([Res[i,3] for i in 1:size(Res,3)])/size(Res,3)
    Res[end,4] = sqrt((sum([Res[i,3] for i in 1:size(Res,3)].-Res[end,3]).^2)/size(Res,3))
    
    path = NamePath(CSVpathNames[1])
    
    CSV.write(path*Title*".csv",Tables.table(Res),header=["AvgKp","StdKp","AvgSize","StdSize"])

    return# AvgKp, StdKp, AvgSize, StdSize
    
end

########################################################################
########################################################################