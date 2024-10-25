function transpiration_I2024(rsd, rld, albedo, tao_s, tao_l, ta, tw, vpd, ra, rc, fq = 1.0)
 
    # This program calculates Et (transpiration mm hr-1)
    # Please refer to Ikawa et al (2024, PCE) for the explanation of the variables.
    # Please feel free to contact hikawa.biomet@gmail.com for any questions.
    # fq: factor for incoming radiation

    sb = 5.67e-8
    md = 28.9
    mv = 18.0
    cp = 1.01e3*md*1e-3 # specific heat of air (Jmol-1K-1)
    lv=2.5e6*mv*1e-3 # latent heat for vaporization Jmol-1 
    T0 = 273.15
    emis_c = 1.0
    emis_w = 1.0

     rsc = (1 .- albedo .- tao_s) .* rsd
     rlc = (1 .- tao_l) .* (rld .+ emis_w * sb * (tw .+ T0) .^4 - 2 * emis_c * sb * (ta .+ T0) .^4)

     dt = ( fq * (rsc .+ rlc) .- lv.*vpd./(press.*(ra .+ rc)) ) ./
    ( cp./ra .+ lv*(deslowe.(ta)*1e-1)./press./(ra .+ rc) .+
    8 * (1 .- tao_l) .* emis_c * sb * (ta .+ T0).^3 )

    Dl = vpd .+ deslowe.(ta)*1e-1 .* dt

    Tc = ta .+ dt
    CeU = 1. ./ (press.*(ra .+ rc))
    LEt = CeU .* Dl * lv # W m-2
    Et =  LEt ./lv * mv * 1e-3 * 3600 # mm hr-1
   
    return LEt, Et, Tc, Dl, CeU

# Checking the reasonability of the ablbedo value based on the canopy property
    # rf = 0.3 # leaf reflectivity
    # tf = 0.2 # leaf transmissivity
    # rg = 0.1
    # df = 1.66
    # lai = 4.5
    # f = 0.5
    # fs = f*df

    # f0=(1-tf)*fs
    # rf0=rf/(1-tf)

    # alpha0=sqrt(1-rf0^2)
    # e=exp(alpha0*f0*lai)
    # ff1=(1-rg*rf0+alpha0)*e
    # ff2=(1-rg*rf0-alpha0)/e
    # a=ff1-ff2
    # b=ff1+ff2
    # albedo=(1-alpha0*b/a)/rf
    # absog=2*alpha0*(1-rg)/a
    # absoc=1-albedo-absog
    # xml=exp(-f*lai*df)

   end

function eslowe(tc)
    #     Lowe, P. R. (1977) polynominal expression.
#     Tc : temperature in (ﾟC)
#     deslow(Tc) = d(es(Tc))/dTc
#     es(Tc) : saturated vapour pressure for 'water' at tc (ﾟC) in (hPa)
#                good aprox. for -50ﾟC .. +70ﾟC (maybe)
#
    es = 6.107799961+tc*(4.436518521*(10^-1)+tc*(1.428945805*(10^-2)
    +tc*(2.650648471*(10^-4)+tc*(3.031240396*(10^-6)+tc*(2.034080948*(10^-8)+tc
    *6.136820929*(10^-11))))))
    return es
end

function deslowe(tc)
   des = 4.438099984*(1e-1)+tc*(2.857002636*(1e-2)+tc*(7.938054040*(1e-4)
    +tc*(1.215215065*(1e-5)+tc*(1.036561403*(1e-7)+tc*(3.532421810*(1e-10)
    -tc*7.090244804*(1e-13))))))
 return des

end
   
function photosynthesis_I2024(co2, vpdl, ta, tl, apar, c_alpha, c_beta, c_gamma, optgs=1)

    # return photosynthesis (A, umolm-2s-1) and stomatal conductance (gsw, molm-2s-1)

 A = zeros(length(ta))
 gsw = zeros(length(ta))

 for i = 1:length(ta)
    (A[i],gsw[i])=
     photosynthesis_I2024_main(apar[i],co2[i],vpdl[i],tl[i],c_alpha, c_beta, c_gamma,optgs)

 end

return A, gsw
end

function photosynthesis_I2024_main(apar,ca,vpdl,tl,c_alpha,c_beta,c_gamma,optgs=1)

# FvCB model with predetermined parameters for plants used in Ikawa et al (2024, PCE)
# FvCB for rice details found in Ikawa et al (2019, PPS) and Ikawa et al (2018, GCB)
# optgs = 1 (Leuning model), 2 (Ball model), 3 (Medlyn model)

    fb = 0.53
    rjv = 1.05    
    gsmin = 0.0001
    vcmax25 = 117. # based on the calculation of (Nl - Nb)*kai = (90 - 3.72)*1.36
    Rd25 = vcmax25*0.006
    gm25 = 0.11 # based on the calculation of (Nl - Nbgm)*kaigm = (90 - 42.4)*0.0023
    gbw = 5.0
    press = 1000.
    
    # initial values [any numbers]  
    A = 0.1
    Ao = A + 0.1 #umolm-2s-1 initial random value
    count = 0
    
    # Model parameter
    EPS = 0.00001
    c_max = 999         # maximum # of iteration 
    
    # Sharkey [2007] T parameterization based on Bernacchi [2002]
    c = zeros(Float64, 7)
    dH = zeros(Float64, 7)
    c[1] = 18.7145; dH[1] = 46.39 # Rd
    c[2] = 26.355; dH[2] = 65.33 # Vcmax
    c[3] = 17.71; dH[3] = 43.9 # Jmax
    c[4] = 35.9774; dH[4] = 80.99 #kc
    c[5] = 12.3772; dH[5] = 23.72 #ko
    c[6] = 11.187; dH[6] = 24.46 # gstar
    
    # Scafaro [2011] temoerature function for rice  
    c[7]=18.29; dH[7]=45.22 # temp function for gm
    
    b = 0.76+0.018 * tl-3.7 * 10^(-4)*tl .^2
    phips2 = 0.352+0.022 * tl-3.4 * 10^(-4)*tl .^2 # quantum yield [umol/mol]
    co = 210000 # oxygen umol mol-1
    
    # Internal calculation for model input()
    tlk = tl+273.15 # K
    pressh = press * 100. # pressure [Pa]
    R = 8.314
    Rd = Rd25*exp(c[1]-dH[1]*1000/R/tlk) # mitochondrial respiration in umolm-2s-1
    gbc = gbw/(1.6^(2/3)) # boundary conductance for CO2 
    
    kc = exp(c[4]-dH[4]*1000/R/tlk)/pressh*10^6           # umol mol-1
    ko = exp(c[5]-dH[5]*1000/R/tlk)*1000/pressh*10^6 
    gstar = exp(c[6]-dH[6]*1000/R/tlk)/pressh*10^6 
    gmtfun = exp(c[7]-dH[7]*1000/R/tlk)
    
    gm = gm25 * gmtfun
    
    vcmax = vcmax25*exp(c[2]-dH[2]*1000/R/tlk)
    jmax25 = vcmax25*rjv
    jmax = jmax25*exp(c[3]-dH[3]*1000/R/tlk)
    j = ((jmax+fb*phips2*apar)-((jmax+fb*phips2*apar)^2-4*b*(jmax*fb*phips2*apar))^(0.5))/2/b
    
    local gsw, cc, Av, Aj
    while abs(Ao - A) > EPS
        count = count + 1
    
        if optgs == 1
        gsw = c_alpha * A /(max(vpdl + c_beta, 0.0001)) / (ca - c_gamma)
        elseif optgs == 2
            rhc = (0.1*eslowe(tl) - vpdl)/(0.1*eslowe(tl))
        gsw = c_alpha * A * rhc / (ca - c_gamma)
        elseif optgs == 3
        gsw =  A * (1 + c_alpha/sqrt(vpdl)) / (ca - c_gamma)
        end
    
        if gsw < gsmin
            gsw = gsmin
        end
        
        gsc = gsw/1.6 # molar diffusivity of CO2/H2O is 1/1.6
       
        cs = ca-A/gbc
        ci = cs-A/gsc
                
        if ci<gstar+A/gm
            ci=gstar+A/gm
        end
    
        Ao = A # when variables are not vector, no need to use copy()
        
        # Ethier et al. (2006)
        prmv_a = -1/gm
        prmv_b = (vcmax-Rd)/gm+ci+kc*(1+co/ko)
        prmv_c = Rd*(ci+kc*(1+co/ko))-vcmax*(ci-gstar)
        Av = (-prmv_b+sqrt(prmv_b^2-4*prmv_a*prmv_c))/prmv_a/2
    
        prmj_a = -1/gm
        prmj_b = (j/4-Rd)/gm+ci+2*gstar
        prmj_c = Rd*(ci+2*gstar)-j/4*(ci-gstar)
        Aj = (-prmj_b+sqrt(prmj_b^2-4*prmj_a*prmj_c))/prmj_a/2
    
        A = min(Av, Aj)
    
        cc = ci - A/gm
    
        if cc < gstar
            cc = gstar
        end
    
        if count> c_max
            break
        end
    end
    
     return A, gsw
    
     end