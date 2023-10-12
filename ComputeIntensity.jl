########################################################################
########################################################################

using TiffImages, Images

########################################################################
########################################################################

function Segment!(coord::Array{Int64},img_pixel,mask_pixel,Gain::Int64,dilute,dense,blur,intThreshold::Float64)
    
    x,y,z = coord
    
    mask_pixel = Float64(mask_pixel)
    
    # remove blur
    if mask_pixel > intThreshold*2/3
        blur[x,y,z] = 1
        dilute[x,y,z] = 0
        dense[x,y,z] = 0
        return 0.0,0.0,0,0,1
        
    # compute dilute intensity
    elseif mask_pixel < intThreshold/3
        dilute[x,y,z] = 1
        dense[x,y,z] = 0
        blur[x,y,z] = 0
        return Float64(img_pixel)*100/Gain,0.0,1,0,0

    # compute dense intensity
    else
        dense[x,y,z] = 1
        dilute[x,y,z] = 0
        blur[x,y,z] = 0
        return 0.0,Float64(img_pixel)*100/Gain,0,1,0
    end
    
    return
    
end

########################################################################
########################################################################

function CheckBlur(nDilute::Int64,nDense::Int64,nBlur::Int64,IntDense::Float64,IntDilute::Float64,x::Int64,y::Int64)
    
    if nBlur>x*y*0.5 || nDilute+nBlur>x*y*0.95 || nDense+nBlur>x*y*0.95
        kp = 0
        intensity = 0
        intensitydiluterelative = 0
    else
        kp = (IntDense/IntDilute)*(nDilute/nDense)
        intensity = IntDense
        intensitydiluterelative = IntDilute/nDilute
    end
    
    return kp, intensity, intensitydiluterelative
end

########################################################################
########################################################################

function ComputeIntensityPlane(img,mask,Gain::Int64,dilute,dense,blur,z::Int64,intThreshold::Float64)
    
    xy = size(img)
    IntDilute = 0
    IntDense = 0
    nDilute = 0
    nDense = 0
    nBlur = 0

    for y in 1:xy[2]
        for x in 1:xy[1]
            intDil,intDen,dil,den,blr = Segment!([x,y,z],img[x,y,z],mask[x,y,z],Gain,dilute,dense,blur,intThreshold)
            IntDilute+=intDil
            IntDense+=intDen
            nDilute+=dil
            nDense+=den
            nBlur+=blr
        end
    end

    # discard planes with too much blur
    kp,intensity,intensitydiluterelative = CheckBlur(nDilute,nDense,nBlur,IntDense,IntDilute,xy[1],xy[2])
    
    return dense, kp, intensity, intensitydiluterelative, nDilute
    
end

########################################################################
########################################################################

function ComputeIntensity(imgName,maskName,Gain::Int64)
    
    # import image and mask
    img = TiffImages.load(imgName)
    mask = TiffImages.load(maskName)
    
    xyz = size(img)
    dilute = deepcopy(img)
    dense = deepcopy(img)
    blur = deepcopy(img)

    Kp = fill(0.0, xyz[3])
    Intensity = fill(0.0, xyz[3])
    IntensityDiluteRelative = fill(0.0, xyz[3])
    NDilute = fill(0, xyz[3])
    
    intThreshold = Float64(maximum(mask))
    
    for z in 1:xyz[3]
        
        dense,kp,intensity,intensitydiluterelative, nDilute = ComputeIntensityPlane(img,mask,Gain,dilute,dense,blur,z,intThreshold)
        Kp[z] = kp
        Intensity[z] = intensity
        IntensityDiluteRelative[z] = intensitydiluterelative
        NDilute[z] = nDilute
    end
    
    return Kp, Intensity, IntensityDiluteRelative, NDilute, dense

end

########################################################################
########################################################################