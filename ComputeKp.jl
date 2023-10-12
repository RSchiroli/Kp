########################################################################
########################################################################

using CSV, Tables, DataFrames, TiffImages, Images

include("ComputeIntensity.jl")
include("Labelling.jl")

########################################################################
########################################################################

function WriteCSV!(name,allKps::Array{Any,1},sizes::Array{Any,1},Label::Array{Any,1},LabelledPlanes::Array{Any,1},Intensity::Array{Any,1},intensityDiluteRel::Array{Float64,1},nDilute::Array{Int64,1})
    
    # parse the name
    path = split(name,"/")
    filename,_ = split(path[end],".")
    csvname = ""
    for p in path[1:end-1]
        csvname *= p*"/"
    end
    csvname *= filename*".csv"
    
    # write in a matrix for CSV
    len = length(allKps)
    result = []
    
    for z in 1:len
        l = length(allKps[z])
        res = fill(0.0,l,7)
        for i in 1:l
            res[i,1] = LabelledPlanes[z]
            res[i,2] = Label[z][i]
            res[i,3] = sizes[z][i]
            res[i,4] = allKps[z][i]
            res[i,5] = Intensity[z][i]
            res[i,6] = intensityDiluteRel[LabelledPlanes[z]]
            res[i,7] = nDilute[LabelledPlanes[z]]
        end
        push!(result,res)
        
        # write CSV file
        if z != 1
            CSV.write(csvname,Tables.table(res),append=true)
        else
            CSV.write(csvname,Tables.table(res),header=["Plane","Label","Size","Kp","Intensity","DiluteIntensity","nDilute"])
        end
        
    end
    
    return
    
end

########################################################################
########################################################################

function CheckPlane(nPlane::Int64,Z::Int64)

    if nPlane+1>Z
        nPlane = Z-1
        @warn "nPlane+1 > number of planes" nPlane+1
    end
    
    return nPlane
    
end

########################################################################
########################################################################

function ComputeKpPlane(img,label::Dict{Any,Any},intensityDiluteRel::Float64,threshold::Int64)

    Kp = []
    int = []
    sizes = []
    lab = []
    for k in keys(label)
        s = size(label[k],1)
        if s>threshold
            intensityDroplet = Float64(sum([img[pixel[1],pixel[2]] for pixel in label[k]]))
            push!(Kp,intensityDroplet/(s*intensityDiluteRel))
            push!(int,intensityDroplet)
            push!(sizes,s)
            push!(lab,k)
        else
            delete!(label,k)
        end
    end
    
    return Kp,int,sizes,lab
    
end

########################################################################
########################################################################

function ComputeKp(imgName,intensity::Array{Float64,1},intensityDiluteRel::Array{Float64,1},nPlane::Int64,label::Array{Dict{Any,Any},1},threshold::Int64)
    
    img = TiffImages.load(imgName)
    maxIntZ = findmax(intensity)[2]
    minZ = maxIntZ-div(nPlane,2)
    maxZ = maxIntZ+div(nPlane,2)
    if minZ < 1
        minZ = 1
    end
    if maxZ > length(intensity)
        maxZ = length(intensity)
    end
    
    avgKp = []
    stdKp = []
    allKps = []
    allIntensities = []
    Sizes = []
    Lab = []
    LabelledPlanes = []
    
    CheckPlane(nPlane,length(intensity))
    
    plane = 1
    for z in minZ:maxZ
        if !isempty(label[plane])
            Kp,int,sizes,lab=ComputeKpPlane(img[:,:,z],label[plane],intensityDiluteRel[z],threshold)
            push!(avgKp,sum(Kp)/length(Kp))
            push!(stdKp,sqrt(sum((Kp.-avgKp[plane]).^2)/length(Kp)))

            # Reorder sizes and Kps
            p = sortperm(lab)
            Kp = Kp[p]
            int = int[p]
            sizes = sizes[p]
            lab = lab[p]

            push!(allKps,Kp)
            push!(allIntensities,int)
            push!(Sizes,sizes)
            push!(Lab,lab)
            push!(LabelledPlanes,z)

            # Relabelling!(label[plane])
            plane += 1
        else
            ################################### problem there
            push!(allKps,[0])
            push!(allIntensities,[0])
            push!(Sizes,[0])
            push!(Lab,[0])
            push!(LabelledPlanes,z)
        end        
    end
    
    return avgKp,stdKp,allKps,Sizes,Lab,allIntensities,LabelledPlanes ### <--- In this version of the code avg and std are returned empty
    
end

########################################################################
########################################################################

function AnalyseZStack(imgName,maskName,nPlane::Int64,Gain::Int64,threshold::Int64)    

    kp, intensity, intensityDiluteRel, nDilute, dense =  ComputeIntensity(imgName,maskName,Gain)

    nPlane = CheckPlane(nPlane,length(kp))

    label = Labelling(dense,intensity,nPlane)

    avgKp, stdKp, allKps, Sizes, Lab, Intensity, LabelledPlanes = ComputeKp(imgName,intensity,intensityDiluteRel,nPlane,label,threshold)

    return avgKp, stdKp, allKps, Sizes, Lab, Intensity, intensityDiluteRel, nDilute, LabelledPlanes, label
        
end
    

########################################################################

function CalcKp(imgName,maskName;nPlane::Int64=0,Gain::Int64=100,threshold::Int64=0,returnAll::Bool=false,returnLabel::Bool=false)
    
    avgKp, stdKp, allKps, Sizes, Lab, Intensity, intensityDiluteRel, nDilute, LabelledPlanes, label = AnalyseZStack(imgName,maskName,nPlane,Gain,threshold)
    
    WriteCSV!(imgName,allKps,Sizes,Lab,LabelledPlanes,Intensity,intensityDiluteRel,nDilute)
    
    (returnAll != false) && return avgKp, stdKp, allKps, Sizes, label, nDilute
    (returnLabel != false) && return avgKp, stdKp, label
    
    return avgKp, stdKp

end

########################################################################
########################################################################
