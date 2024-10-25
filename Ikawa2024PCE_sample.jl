using CSV
using DataFrames
using Plots
include("C:.\\Ikawa2024PCE_subroutine.jl")

# This program calculates canopy transpiration based on Ikawa et al (2024, PCE)

# Et: canopy transpiration (mm hr-1)
# Tc: canopy temperature (deg C)
# Dl: canopy - air vapor pressure deficit (kPa)

global Et, Tc, Dl

sb = 5.67e-8
md = 28.9
mv = 18.0

albedo = 0.2
tao_s = 0.05
tao_l = 0.02 

cp = 1.01e3*md*1e-3 # specific heat of air (Jmol-1K-1)
lv=2.5e6*mv*1e-3 #  latent heat for vaporization Jmol-1 
T0 = 273.15
press = 100.

ga_high = 7.6 # mol m-2 s-1
ga_low = 3.3 # mol m-2 s-1

c_alpha1 = 3.73 # The greater the smaller Rc
c_alpha2 = 16.2
c_beta = 1.18
c_gamma = 106.

# read daily data from the TGC experiment
fid = ".\\forcing\\TGCsummary_hourly.csv"

data=CSV.File(fid,header=1,skipto=2,delim=","; missingstrings=["-9999", "NA"]) |> DataFrame

# variables to define the conditions
doy = data.doy
hour = data.hour
Treatment = data.Treatment

# model forcing variables
rsd = data.SR # downward shortwave radiation (W m-2)
rld = data.Ld # downward longwave radiation (W m-2)
ta = data.Ta # air temperature (deg C)
vpd = data.vpd # vapor pressure deficit (kPa)
co2 = data.CO2data # CO2 concentration (umol mol-1)
tw = data.Tw # water surface tempreature (deg C)

# observation variables for the validations
Etobs = data.Et # canopy transpiration
Tcobs = data.Tc # canopy temperature

rsc = (1 .- albedo .- tao_s) .* rsd # absorbed shortwave radiation
apar = rsc * 2.0 # absorbed par

esat = eslowe.(ta) # hpa
ea = esat .- vpd*10. # hpa

fnddata_Ctrl = (((Treatment .== "Control") .| (Treatment .== "e-T")) .& (doy .!= 224) )
fnddata_eC = (((Treatment .== "e-C") .| (Treatment .== "e-TC")) .& (doy .!= 224) )

ga = zeros(length(ta))
ga .= ga_high
ga[(hour .<= 9) .| (hour .== 16)] .= ga_low

ra = 1 ./(ga)

LEtobs = Etobs * lv / mv / 1e-3 / 3600 # mm hr-1 to W m-2
Dlobs = 0.1 * (eslowe.(Tcobs) - ea) 
# Dlobs = vpd .+ deslowe.(ta)*1e-1 .* (Tcobs - ta) # first order approximation in Ikawa et al (2024, PCE)
rcobs = lv .* Dlobs ./ press ./ LEtobs - ra

# initial values
Dl = copy(vpd)
Tc = copy(ta)

for l = 1:2 # minimum iteration to utilize calculated Tc and Dl in the photosynthesis model
~, gsw = photosynthesis_I2024(co2, Dl, ta, Tc, apar, c_alpha2, c_beta, c_gamma)
rc = 1. ./ (gsw*c_alpha1)
~, Et, Tc, Dl = transpiration_I2024(rsd, rld, albedo, tao_s, tao_l, ta, tw, vpd, ra, rc)
end

scatter(Etobs, Et)