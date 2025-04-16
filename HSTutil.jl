#!/usr/bin/env julia
using FITSIO, Plots, Printf, NumericalIntegration, StatsBase, DelimitedFiles, DSP, Interpolations, CubicSplines, Distributions
#include("pyimports.jl")
# using .pyimports: numpy
#sorting all visit files
#files = readdir("/home/kirk/Documents/research/Dexter/STORM/HST/models/")


################################# DATA STRUCTURES AND READING FUNCTIONS ################################
function getSortedFiles(spectra_dir=pwd())
    files = readdir(spectra_dir)
    function file_sort_key(file_name)
        sub_string = split(file_name, "_")[5]
        file_num = split(sub_string, "-")[end] #visit "number"
        if isnumeric(file_num[3])
            # If the file name ends with a number, sort numerically
            return parse(Int, file_num[2:end])
        else
            # If the file name ends with a letter, add letter position in alphabet to get numerical position (numbers stop at 99)
            return (parse(Int, file_num[2])) + 100 + 26 * (parse(Int, file_num[2])) + Int(file_num[end]) - Int('A') + 1
        end
    end
    sorted_files = sort(files, by=file_sort_key)
    return sorted_files
end

function getSpectrumTimeInfo(file)
    f = FITS(file)
    header_string = read_header(f[1],String)
    exp = split(split(header_string,"EXPTIME =")[2])[1] #s
    start = split(split(header_string,"HJDSTART=")[2])[1] #heliocentric julian date
    return parse(Float64,exp), parse(Float64,start)
end

struct spectrum
    filename::String
    λ::Array{Float64,1}
    flux::Array{Float64,1}
    error::Array{Float64,1}
    cont_flux::Union{Array{Float64,1},Nothing}
    model_flux::Union{Array{Float64,1},Nothing}
    t::Float64
    exp::Float64
end

struct lcDatum #one piece of light curve
    t::Float64 #heliocentric julian Date
    F1158::Float64 #flux at 1158 A
    e_F1158::Float64 #1 sigma error in flux at 1158 A
    F1367::Float64 #flux at 1367 A
    e_F1367::Float64 #1 sigma error in flux at 1367 A
    F1469::Float64 #flux at 1469 A
    e_F1469::Float64 #1 sigma error in flux at 1469 A
    F1745::Float64 #flux at 1745 A
    e_F1745::Float64 #1 sigma error in flux at 1745 A
    F_Lya::Float64 #flux of broad Ly alpha emission
    e_F_Lya::Float64 #1 sigma error in Ly alpha emission
    F_NV::Float64 #flux of broad NV emission
    e_F_NV::Float64 #1 sigma error in NV emission
    F_SiIV::Float64 #flux of broad SiIV emission
    e_F_SiIV::Float64 #1 sigma error in SiIV emission
    F_CIV::Float64 #flux of broad CIV emission
    e_F_CIV::Float64 #1 sigma error in CIV emission
    F_HeII::Float64 #flux of broad HeII emission
    e_F_HeII::Float64 #1 sigma error in HeII emission
end

struct LC
    t::Array{Float64,1} #days
    flux::Array{Float64,1} #flux
    error::Array{Float64,1} #1 sigma error in flux
end

struct LC2D
    LCs::Array{LC,1}
    λ::Array{Float64,1}
end

function getLCData(file)
    fileData = readdlm(file,skipstart=30) #first 30 lines are header info pasted above
    data = Array{lcDatum,1}(undef, size(fileData,1)) #initialize data Array
    for i = 1:size(fileData,1)
        data[i] = lcDatum(fileData[i,1], fileData[i,2], fileData[i,3], fileData[i,4], fileData[i,5], fileData[i,6], fileData[i,7], fileData[i,8], fileData[i,9], fileData[i,10], fileData[i,11], fileData[i,12], fileData[i,13], fileData[i,14], fileData[i,15], fileData[i,16], fileData[i,17], fileData[i,18], fileData[i,19])
    end
    return data
end

function getLC(LCData,name,tStop=75.)
    t = [LCData[i].t for i in 1:length(LCData)]
    t = t.-t[1] #normalize to first visit
    if name == "Ly alpha"
        name = "Lya"
    end
    if tryparse(Int,name) == nothing
        name = "F_"*name
    else
        name = "F"*name
    end
    f = [getproperty(LCData[i],Symbol(name)) for i in 1:length(LCData)]
    e = [getproperty(LCData[i],Symbol("e_"*name)) for i in 1:length(LCData)]
    mask = t.<tStop
    return LC(t[mask],f[mask],e[mask])
end

function getSpectra(sortedFiles; src_dir=pwd(),model="all")
    data = Array{spectrum,1}(undef, length(sortedFiles))
    for (i,file) in enumerate(sortedFiles)
        filename = string(src_dir,file)
        f = FITS(filename)
        λ = model != false ? read(f[2], "WAVE") : read(f[2], "WAVELENGTH") #angstroms
        flux = read(f[2], "FLUX") #cW / m^2 / A
        error = model != false ? read(f[2], "e_FLUX") : read(f[2], "ERROR") #"          "
        cont_flux = model != false ? read(f[2], "FCONT") : nothing  #"      "
        model_flux = nothing
        if model != false
            if model == "all"
                model_flux = read(f[2], "MODFLUX")
            else
                if typeof(model) == Array{Int64,1}
                    model_flux = read(f[2],"GAUSS$(model[1])")
                    for j in model[2:end]
                        model_flux .+= read(f[2],"GAUSS$j")
                    end
                elseif typeof(model) == Int64 || typeof(model) == Int32
                    model_flux = read(f[2],"GAUSS$model")
                else
                    error("model must be an integer or array of integers")
                end
            end
        end
        exp,t = getSpectrumTimeInfo(filename)
        data[i] = spectrum(filename, λ, flux, error, cont_flux, model_flux, t, exp)
    end
    return data
end

function trapint(x,y)
    int = 0.0
    for i in 1:length(x)-1
        int += (y[i]+y[i+1])*(x[i+1]-x[i])/2
    end
    return int
end

function get2DLC(vData=nothing;LCfile="hlsp_storm_hst_cos_ngc-5548_g130m-g160m_v1_emission-lines.txt", spectra_dir="HST/models/", line="Ly alpha", cont="1367", tStop=75., fig5File="HST/Lya_T1_VBinnedDelay_FromFigure5.csv",λCen=1235.,useModel=true)
    sortedfiles = getSortedFiles(spectra_dir)
    models = getLineModelStrings(line)
    lineData = getSpectra(sortedfiles,src_dir=spectra_dir,model=models)
    λUnique = sort(unique(vcat([spectrum.λ for spectrum in lineData]...)))
    vUnique = (λUnique .- λCen)./λCen .* 3e5
    if isnothing(vData)
        vData = vUnique
    end
    diffs = [vData[i+1]-vData[i] for i in 1:length(vData)-1]
    diff = mean(diffs)
    vEdges = [vData[1]-diff/2;vData.+diff/2]
    LCs = Array{LC}(undef,length(vData))
    ts = [spectrum.t for spectrum in lineData]
    for vi in eachindex(vData)
        fluxes = zeros(length(vData),length(ts))
        errors = zeros(length(vData),length(ts))
        for (ti,spectrum) in enumerate(lineData)
            v = (spectrum.λ .- λCen)./λCen .* 3e5
            mask = (v .> vEdges[vi]) .& (v .< vEdges[vi+1])
            if sum(mask) > 0
                flux = useModel ? spectrum.model_flux[mask] : spectrum.flux[mask]
                error = spectrum.error[mask]
                t = spectrum.t
                fluxes[vi,ti] = trapint(spectrum.λ[mask],flux)
                errors[vi,ti] = sqrt(trapint(spectrum.λ[mask],error.^2))
            end
        end
        LCs[vi] = LC(ts,fluxes[vi,:],errors[vi,:])
    end
    λs = vData./3e5 .* λCen .+ λCen
    return LC2D(LCs,λs)
end

struct avgSpectrum
    λ::Array{Float64,1}
    flux::Array{Float64,1}
    error::Array{Float64,1}
    model_flux::Array{Float64,1}
end

function getAvgSpectrum(spectra; λRange=nothing, tStop=75.)
    avg_λ = sort(unique(vcat([spectra[i].λ for i in 1:length(spectra)]...)))
    if λRange != nothing
        mask = (avg_λ .> λRange[1]) .& (avg_λ .< λRange[2])
        avg_λ = avg_λ[mask]
    end
    avg_flux = zeros(length(avg_λ))
    avg_error = zeros(length(avg_λ))
    avg_model_flux = zeros(length(avg_λ))
    n = zeros(length(avg_λ))
    t0 = minimum([spectrum.t for spectrum in spectra])
    for (i,λ) in enumerate(avg_λ)
        for spectrum in spectra
            if (λ in spectrum.λ) && ((spectrum.t-t0) < tStop)
                ind = findfirst(spectrum.λ .== λ)
                avg_flux[i] += spectrum.flux[ind]
                avg_error[i] += spectrum.error[ind]^2
                avg_model_flux[i] += spectrum.model_flux[ind]
                n[i] += 1
            end
        end
    end
    return avgSpectrum(avg_λ,avg_flux./n,sqrt.(avg_error)./n,avg_model_flux./n)
end

function getVariability(spectra,avgSpectrum; tStop=75.)
    avg_λ = avgSpectrum.λ
    variability = zeros(length(avg_λ))
    n = zeros(length(avg_λ))
    t0 = minimum([spectrum.t for spectrum in spectra])
    for (i,λ) in enumerate(avg_λ)
        for spectrum in spectra
            if (λ in spectrum.λ) && ((spectrum.t-t0) < tStop)
                ind = findfirst(spectrum.λ .== λ)
                variability[i] += (spectrum.flux[ind] - avgSpectrum.flux[i])^2
                n[i] += 1
            end
        end
    end
    return sqrt.(variability./n)
end

function getLineModelStrings(line::String; sample_fits::String="/home/kirk/Documents/research/Dexter/STORM/HST/models/hlsp_storm_hst_cos_ngc-5548-go13330-v0a_g130m-g160m_v1_model.fits") #don't think this is quite right at end product, do further testing tmoorrow
    f = FITS(sample_fits)
    header_string = read_header(f[1],String)
    split_strip_header = strip.(split(header_string,"HISTORY"))
    GaussianModelStrings = split_strip_header[startswith.(split_strip_header,"Gaussian")]
    lineMask = occursin.(line,GaussianModelStrings)
    lineStrings = GaussianModelStrings[lineMask]
    BLRMask = occursin.("Broad",lineStrings)
    BLAMask = occursin.("Broad Absorption",lineStrings)
    BLREmissionMask = (BLRMask) .& (.!BLAMask)
    bumpMask = occursin.("Bump",lineStrings)
    finalMask = (BLREmissionMask) .|| (bumpMask)
    lineStrings = lineStrings[finalMask]
    components = [strip(split(lineStrings[i]," ")[3],':') for i in eachindex(lineStrings)]
    return tryparse.(Int,components)
end

function getHSTData(LCfile,spectra_dir,line="CIV",cont="1367",tStop=75.;wrap=true)
    λRange = nothing
    if line == "CIV"
        λRange = [1450,1725]
    elseif line == "Ly alpha"
        λRange = [1190,1275]
    else
        error("line not recognized or λRange not specified -- see getHSTData for default values")
    end
    LCData = getLCData(LCfile)
    lineLC = getLC(LCData,line,tStop)
    CLC = getLC(LCData,cont,tStop)
    sortedfiles = getSortedFiles(spectra_dir)
    models = getLineModelStrings(line)
    lineData = getSpectra(sortedfiles,src_dir=spectra_dir,model=models)
    avgSpectrum = getAvgSpectrum(lineData,λRange=λRange,tStop=tStop)
    if wrap
        return lineLC.t,lineLC.flux,lineLC.error,CLC.t,CLC.flux,avgSpectrum.λ,avgSpectrum.model_flux,avgSpectrum.error
    else
        return lineLC,CLC,avgSpectrum
    end
end

# function getHSTDataWrap(input)
#     LCfile,spectra_dir,line,cont,tStop = input
#     if typeof(tStop) != Float64
#         tStop = tryparse(Float64,tStop)
#     end
#     line_t,lineLC,LCerror,cont_t,contLC,λ,model_flux,LPerror = getHSTData(LCfile,spectra_dir,line,cont,tStop)
#     return numpy[].array([line_t,lineLC,LCerror,cont_t,contLC,λ,model_flux,LPerror],dtype="object")
# end

############################# ANALYSIS FUNCTIONS / PLOTTING #############################
#do polynomial interpolation of LC data to fill in between timestamps, specify N points between each tBinned
function getSpline(x::Array{Float64,1},y::Array{Float64,1})
    if length(x) != length(y)
        error("x and y must be same length")
        return nothing
    end
    spline = CubicSpline(x,y)
    f(x) = spline[x]
    return f
end


function plotLC(LCData, lineList; tRange = nothing, spline = false, model = nothing)
    t = [LCData[i].t for i in 1:length(LCData)]
    t = t.-t[1] #normalize to first visit
    mask = [true for i in 1:length(LCData)]
    p = plot(title = "HST STORM Light Curves",xlabel="t [days after first visit]", ylabel="Flux [normalized (1e-15 continuum; 1e-13 line)]", legend=:topleft)
    if tRange != nothing
        mask = (t.>tRange[1]) .& (t.<tRange[2])
        t = t[mask]
    end
    for (i,line) in enumerate(lineList)
        f = [getproperty(LCData[i],Symbol(line)) for i in 1:length(LCData)]
        norm = line == "F1367" ? 1e-15 : 1e-13 #normalize to 1e-15 for continuum, 1e-13 for emission lines
        e = [getproperty(LCData[i],Symbol("e_"*line)) for i in 1:length(LCData)]
        p = plot!(t,f[mask]./norm,ribbon=e[mask]./norm,marker=:circle,markerstrokewidth=0.,label=line,lw = spline ? 0. : 2.,linealpha = spline ? 0. : 1.,c=i)
        if spline
            interp = getSpline(t,f[mask]./norm)
            tInterp = range(t[1],t[end],length=1000)
            p = plot!(tInterp,interp.(tInterp),label="",lw=2.,linealpha=1.,c=i)
        end
        if model != nothing && line != "F1367" #assume model = [model_continuum_LC, model_emission_LC]
            data_span = maximum(f) - minimum(f)
            f_model = (model[2].*data_span .+ f[1])./norm #zero point is first flux value, match span (initial span should be 1)
            if length(t) != length(f_model)
                maxInd = minimum([length(t),length(f_model)])
                t_model = t[1:maxInd]; f_model = f_model[1:maxInd]
            end
            p = plot!(t_model,f_model,marker=:star,markerstrokewidth=0.,label=line*" model",lw = spline ? 0. : 2.,linealpha = spline ? 0. : 1.,c=i,ls=:dash)
            if spline
                interp = getSpline(t_model,f_model)
                tInterp = range(t_model[1],t_model[end],length=1000)
                p = plot!(tInterp,interp.(tInterp),label="",lw=2.,linealpha=1.,c=i,ls=:dash)
            end
        end
    end
    p = plot!(ylims=(0,100),size=(720,720),minorticks=true,widen=false)
    return p
end

function syntheticLC(Ct,CLC,Ψτ;continuum="1367",tStop=75.,spline=false,tStart=0.0) #assumes Ψτ is already interpolated to continuum timestamps as in functions.jl/getProfiles
    #Ψτ = Ψτ./maximum(Ψτ) #normalize
    #LC(t) = integral_0^t Ψτ(t-τ)*C(τ)dτ where C(τ) is the continuum light curve
    #write function that does integral above and use to get synthetic light curve
    #C = getLC(LCData,continuum)
    #C needs to be interpolated to match tBinned -- actually instead let's interpolate Ψτ and tBinned to match continuum t
    #in reality later in modelling we can just sample tBinned at continuum t but for simple exercise do this
    t = Ct
    mask = t.<tStop
    t = t[mask]; C = CLC[mask]; Ψτ = Ψτ[mask] #only go to tStop
    C = C./maximum(C)
    span = maximum(C) - minimum(C) 
    ΔC = C .- C[1]
    ΔC = ΔC ./ span #make span 1, normalize to first point so that + = brighter and - = dimmer
    if spline 
        CSpline = getSpline(t,ΔC)
        ΨτSpline = getSpline(t,Ψτ)
        t = range(t[1],t[end],length=1000)
        ΔC = CSpline.(t); Ψτ = ΨτSpline.(t)
    end
    function LC(t,Ψτ,ΔC)
        N = length(t)
        LC = zeros(N)
        for ti=1:N
            integrand = 0
            for τ = 1:ti
                dτ = ti < N ? t[τ+1]-t[τ] : t[end] - t[end-1]
                integrand += Ψτ[τ]*ΔC[ti-τ+1]*dτ
            end
            LC[ti] = integrand
        end
        return LC #this is really like ΔLC
    end
    Δlc = LC(t,Ψτ,ΔC) #these are variations, not properly normalized
    #zeroPoint = findfirst(t.>tStart) #find first index after tStart, this is value that should be "zero" 
    #Δlc = Δlc .- Δlc[zeroPoint] #normalize to zero point
    #spanlc = maximum(Δlc) - minimum(Δlc) #when normalized the extent of variations in continuum should be the same as in emission line
    return Δlc#./spanlc,ΔC #checked by doing CCF and this produces same thing as data!
end

#now generate cross correlation function from getSpectra
function getCCF(u,v,tu,tv; z = 0.0, tRange = nothing, normalize=true, lags = nothing, superSample = 100)
    #u and v are LC to cross correlate (i.e. u = CLC and v = lineLC), t is shared time array
    #lags can either be an integer array of index shifts to try, a minimum and maximum time shift, or nothing (default)
    #returns an array of cross correlation values
    #note that this assumes u and v are sampled at the same times in t, see note below

    #TO-DO -- add interpolation in case spectra t and continuum t don't match up -- or is this already done? confused by paper vs. txt file
    umask = Bool.(ones(Int,length(tu))); vmask = Bool.(ones(Int,length(tv)))
    tu = tu.-tu[1]; tv = tv.-tv[1] #normalize to first visit
    tu = tu./(1+z); tv = tv./(1+z) #rest frame of emitter, correct for redshift
    tList = [tu,tv]
    if length(tu) == length(tv)
        if sum(tu.==tv) == length(tv)
            tList = [tu]
        end
    end
    if length(tu) != length(u) || length(tv) != length(v)
        error("tu/tv and u/v must be same length")
        return nothing
    end
    if tRange != nothing
        umask = (tu.>tRange[1]) .& (tu.<tRange[2])
        vmask = (tv.>tRange[1]) .& (tv.<tRange[2])
        tu = tu[umask]; tv = tv[vmask]
    end
    CCFList = Array{Array{Float64,1},1}(undef,length(tList)); lagsList = Array{Array{Float64,1},1}(undef,length(tList))
    ui = copy(u); vi = copy(v)
    for (i,t) in enumerate(tList)
        tLags = range(t[1],t[end],length=superSample*length(t)) #super sample t to get more accurate interpolation
        uInterp = LinearInterpolation(tu,u[umask],extrapolation_bc=Line())
        vInterp = LinearInterpolation(tv,v[vmask],extrapolation_bc=Line())
        ui = uInterp.(tLags); vi = vInterp.(tLags) #interpolate to super sampled t, so spacing in lags is fixed
        if normalize
            ui = ui./maximum(ui)
            vi = vi./maximum(vi)
        end
        if lags == nothing || (typeof(lags) != Array{Float64,1} && typeof(lags) != Array{Int64,1}) #default lags from docs
            lags = collect(-minimum(Int,(size(ui,1)-1, floor(10*log10(size(ui,1))))):minimum(Int,(size(ui,1), floor(10*log10(size(ui,1))))))
        elseif length(lags) == 2 #lagRange
            minInd = -findfirst(tLags.>abs(lags[1]))
            maxInd = findfirst(tLags.>lags[2])
            lags = collect(minInd:maxInd)
        end
        CCF = crosscor(ui,vi,lags) #cross correlation function
        CCFList[i] = CCF
        lagsList[i] = vcat(-reverse(tLags[1:abs(lags[1])+1]),tLags[2:lags[end]+1])
    end
    #lags[1] = 0 -> t[1]; lags[end] = max +shift -> t[lags[end]+1]; lags[1] = max -shift -> -t[lags[1]+1]
    return mean(CCFList),mean(lagsList) #lags here are index offsets from 0
end #trying to reproduce figure 4 in https://iopscience.iop.org/article/10.1088/0004-637X/806/1/128/pdf , think it works now

getCCFPeak(CCF,lags) = lags[findmax(CCF)[2]] #peak of CCF
function getCCFCentroid(CCF,lags;rCut = 0.8)
    rMax = maximum(CCF)
    mask = CCF.>rCut*rMax
    return  sum(CCF[mask].*lags[mask])/sum(CCF[mask]) #centroid of CCF with cutoff imposed
end

function MCLC(LC,err,t) #get Monte Carlo light curve following FR/FSS outlined in appendix here: https://iopscience.iop.org/article/10.1086/423269/pdf
    N = length(LC)
    d = DiscreteUniform(1,N)
    inds = rand(d,N)
    uniqueCounts = countmap(inds)
    uniqueInds = unique(inds)
    MCLC = zeros(length(uniqueInds))
    tMC = t[uniqueInds]
    for (MCLCi,LCi) in enumerate(uniqueInds)
        d = Normal(LC[LCi],err[LCi]/sqrt(uniqueCounts[LCi]))
        MCLC[MCLCi] = rand(d)
    end
    perm = sortperm(tMC) #put the points back in ascending order in t
    return MCLC[perm], tMC[perm]
end

function getDataLC(lineData,Δλ; λRange=[1520,1646],tStop = nothing, contRange = nothing, model=true)
    # THIS ISNT WORKING RIGHT -- need to remove continuum? see note below, maybe just re-write integration routine seems too complex as is
    #TO DO: "de-trend" and remove continuum level from line (use points at either side of λRange?)
    tLC = [lineData[i].t for i in 1:length(lineData)]
    tLC = tLC.-tLC[1] #normalize to first visit
    if tStop != nothing
        stopInd = findfirst(tLC.>=tStop)
        lineData = lineData[1:stopInd]
        tLC = tLC[1:stopInd]
    end
    if contRange == nothing
        contRange = [(λRange[1]-10,λRange[1]), (λRange[2],λRange[2]+10)]
    end

    function integrateSpectrum(spectrum,λRange,contRange,model)
        f = spectrum.flux; λ = spectrum.λ; error = spectrum.error
        if model
            f = spectrum.model_flux
        end
        λ1,λ2 = mean.(contRange)
        f1,f2 = mean(f[(λ.>contRange[1][1]) .& (λ.<contRange[1][2])]),mean(f[(λ.>contRange[2][1]) .& (λ.<contRange[2][2])])
        m = (f2-f1)/(λ2-λ1); b = f1 - m*λ1
        C(λ,m,b) = m*λ + b
        int = 0.0; cumErr = 0.0
        integrand(λ1,λ2,f1,f2,m,b) = (λ2-λ1)*((f1-C(λ1,m,b))+(f2-C(λ2,m,b)))/2 #trapezoidal rule
        mask = (λ.>λRange[1]) .& (λ.<=λRange[2])
        λ = λ[mask]; f = f[mask]
        for λi in 2:length(λ)
            int += integrand(λ[λi-1],λ[λi],f[λi-1],f[λi],m,b)
            cumErr += ((error[λi]+error[λi-1])*(λ[λi]-λ[λi-1])/2)^2
        end
        return int, sqrt(cumErr)
    end

    nBins = Int(ceil((λRange[2]-λRange[1])/Δλ))
    λBinEdges = range(λRange[1],λRange[2],length=nBins+1) 
    LC = zeros(length(tLC),nBins); LCerror = zeros(length(tLC),nBins)

    for (ti,spectrum) in enumerate(lineData)
        for λi in 1:length(λBinEdges)-1
            λRange = (λBinEdges[λi],λBinEdges[λi+1])
            int, err = integrateSpectrum(spectrum,λRange,contRange,model)
            LC[ti,λi] = int; LCerror[ti,λi] = err
        end
    end
    return LC,LCerror,λBinEdges,tLC
end

#scratch space
function readFigureData(csv)
    data = readdlm(csv,',') #shape [nPoints,5] -- [v,delay,pointNumber,identifier (i.e. Std),second identifier (i.e. +/- 1)]
    v = zeros(Int(size(data)[1]/3)); delay = zeros(Int(size(data)[1]/3)); errUp = zeros(Int(size(data)[1]/3)); errDown = zeros(Int(size(data)[1]/3))
    for i=1:length(v)
        v[i] = mean(data[3*i-2:3*i,1])
        delay[i] = data[3*i-1,2]
        errDown[i] = data[3*i-2,2] - delay[i] #note that then errDown < 0, errUp > 0
        errUp[i] = data[3*i,2] - delay[i] #when plotting then need to take abs(errDown) for yerr/ribbon args 
    end
    return v,delay,errUp,errDown
end

function getLineDelays(lineData,CLC;Δλ=10.,λRange=[1450,1730],model=true,centroid=true,MCMC=1,kwargs...)
    tC = CLC.t; u = CLC.flux;
    LC,LCerr,λBins,tLC = getDataLC(lineData,Δλ,tC,λRange=λRange,model=model)
    lineDelays = zeros(size(LC,2),MCMC)#one delay for each bin
    for MCi in 1:MCMC
        percentFloor = floor(Int,MCi/MCMC*100)
        progressBar = percentFloor > 0 ? "["*"="^(percentFloor-1)*">"*" "^(99-(percentFloor-1))*"]" : "["*">"*" "^(99)*"]"
        print("MC estimation progress: "*progressBar*"\r")
        for i in 1:size(LC,2)
            v = LC[:,i]
            MCLCi,tLCi = MCLC(v,LCerr[:,i],tLC)
            CCF, lags = getCCF(u,MCLCi,tC,tLCi; kwargs...)
            lineDelays[i,MCi] = centroid ? getCCFCentroid(CCF,lags) : getCCFPeak(CCF,lags)
        end
    end
    println()
    binCenters = λBins[1:end-1] .+ Δλ/2
    lineDelays = mean(lineDelays,dims=2)
    return lineDelays, binCenters
end

function plotCCF(u,v,t; tRange=nothing,normalize=true,lags=nothing,spline=false,superSample=100)
    CCF, t = getCCF(u,v,t,tRange=tRange,normalize=normalize,lags=lags,spline=spline,superSample=superSample)
    p = plot()
    if length(CCF) == 1
        p = plot!(t,CCF,legend=false,xlabel="offset τ [days]",ylabel="r",ylims=(-0.1,1),widen=false,size=(720,720),minorticks=true,
            linewidth=2.,c=3,marker=:circle,markerstrokewidth=0.,title="$emission delays relative to $continuum",label="",yticks=[0.2*i for i=0:5])
        τpeak = t[findmax(CCF)[2]]
        p = vline!([τpeak],label="peak delay = $(round(τpeak,digits=2)) days",c=:crimson,lw=2,ls=:dash,legend=:topright)
    else
        labels=["data","model"]
        for (i,ccf) in enumerate(CCF)
            p = plot!(t,ccf,legend=false,xlabel="offset τ [days]",ylabel="r",ylims=(-0.1,1),widen=false,size=(720,720),minorticks=true,
                linewidth=2.,c=i+2,marker=:circle,markerstrokewidth=0.,title="$emission delays relative to $continuum",label=labels[i],yticks=[0.2*i for i=0:5])
            τpeak = t[findmax(ccf)[2]]
            p = vline!([τpeak],label="peak delay = $(round(τpeak,digits=2)) days",c=i+2,lw=2,ls=:dash,legend=:topright)
        end
    end
    if lags != nothing
        p = plot!(xlims=lags)
    end

    return p
end

function LC_CCF_stacked(LCData; continuum="F1367",emission="F_Lya",tRange=nothing,normalize=true,lags=nothing,spline=false,model=nothing)
    p1 = plotLC(LCData, [continuum,emission],tRange=tRange,spline=spline,model=model)
    p2 = plotCCF(LCData,continuum=continuum,emission=emission,tRange=tRange,normalize=normalize,lags=lags,spline=spline,model=model)
    p = plot(p1,p2,layout=@layout([a;b]),size=(720,1440),left_margin=10*Plots.Measures.mm)
    return p
end
