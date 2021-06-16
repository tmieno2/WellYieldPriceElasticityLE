ASCE_PenmanMonteith <- function(Tmax,Tmin,SolRad,Vp,Wspd,Precip,Year,Month,Day,Latitude,Elevation){
  # ---------------------------------------------------------------------- #
  # Function to calculate daily reference ET0 using the ASCE standardized  #
  # reference evapotranspiration equation (Walter et al., 2005)            #
  # ---------------------------------------------------------------------- #
  # Key references:
  #   1. Allen, R.G. et al. (1998) Crop Evapotranspiration: Guidlines for 
  #      computing crop water requirements. Irrig. and Drain. Paper 56, Food 
  #      and Agriculture Organization of the United Nations, Rome.
  #   2. Hargreaves, G.H. and Samani, Z.A. (1982) Estimating potential
  #      evapotranspiration. J. Irrig. and Drain. Engrg., ASCE, 
  #      108(3):225:230
  #   3. Walter, I.A. et al (2005) The ASCE Standardized Reference
  #      Evapotranspiration Equation.
  #   4. Hargreaves, G.H. (1994) Simplified coefficients for estimating 
  #      monthly solar radiation in North America and Europe.
  
  # R version of Tim Foster's Matlab script ASCE_PenmanMonteith.m
  
  # Define constants
  
  #Atmospheric pressure (kPa)
  AtmP <- 101.3*((293-0.0065*Elevation)/293)^5.26
  
  # Psychometric constant (kPa/C)
  psy <- 0.665*10^(-3)*AtmP  # Added parentheses around -3 - James

  # Albedo for grass reference crop (Allen et al. 1998)
  albedo <- 0.23

  # Assume soil heat flux is negligible
  G <- 0
  
  # Convert latitude from degrees to radians
  LatRad <- (pi/180)*Latitude # Rounding error?

  # Stefan-Boltzmann constant (MJ/K^4/m2/day)
  sbc <- 4.903*10^(-9)  # Added parentheses around -9 - James

  # Solar constant (MJ/m2/min)
  Gs <- 0.0820
  
  # Numerator constant (Walter et al. 2005)
  Cn <- 900 # Default for daily calculations for short-reference grass crop

  # Denominator constant (Walter et al. 2005)
  Cd <- 0.34 # Default for daily calculations for short-reference grass crop
  
  # NOTE - for alfalfa reference crop, Cn and Cd should be set equal to 1600
  # and 0.38 respectively.
  
  
  # Calculate additional weather variables
  # Mean daily temperature (degC)
  Tmean <- (Tmax+Tmin)/2

  # Saturation vapour pressure at maximum temperature
  e0max <- 0.6108*exp((17.27*Tmax)/(Tmax+237.3))

  # Saturation vapour pressure at minimum temperature
  e0min <- 0.6108*exp((17.27*Tmin)/(Tmin+237.3))

  # Saturation vapour pressure
  es <- (e0max+e0min)/2

  # Actual vapour pressure 
  ea <- Vp
  
  # Slope of the saturation vapour pressure curve
  svp_slope <- (4098*(0.6108*exp((17.27*Tmean)/(Tmean+237.3))))/((Tmean+237)^2)
  
  # Calculate Julian days
  J <- Day-32+floor(275*(Month/9))+(2*floor(3/(Month+1)))+floor((Month/100)-(((Year%%4)/4))+0.975)  # Replace mod(Year,4) with %% - James

  # Calculate Penman-Monteith variables 
  # Inverse relative distance Earth-Sun
  dr <- 1+0.033*cos(((2*pi)/365)*J)
  
  # Solar decimation (rad)
  sd <- 0.409*sin(((2*pi)/365)*J-1.39)
  
  # Sunset hour angle (rad)
  ws <- acos((-tan(LatRad))*tan(sd)) # ????
  
  # Extraterrestrial radiation (MJ/m^2/day)
  Ra <- (24*60)/pi*Gs*dr*(ws*sin(LatRad)*sin(sd)+cos(LatRad)*cos(sd)*sin(ws))
  
  # Calculate shortwave radiation if not provided as an input
  Rs <- SolRad
  
  # Net shortwave radiation (MJ/m^2/day)
  Rns <- (1-albedo)*Rs
  
  # Clear-sky solar radiation (MJ/m^2/day)
  Rs0 <- (0.75+(0.00002*Elevation))*Ra
  
  # Cloudiness function
  fcd <- 1.35*(Rs/Rs0)-0.35
  fcd <- ifelse(fcd>1.0,1.0,fcd)
  fcd <- ifelse(fcd<0.05,0.05,fcd)
  
  # Net longwave radiation (MJ/m^2/day)
  Rnl <- sbc*fcd*(0.34-0.14*sqrt(ea))*(((Tmax+273.16)^4+(Tmin+273.16)^4)/2)
  
  # Net radiation (MJ/m^2/day)
  Rn <- Rns-Rnl
  
  # Adjust wind speed measurement
  Wspd <- Wspd
  
  
  # Calculate reference evapotranspiration
  Et0 <- ((0.408*svp_slope*(Rn-G))+psy*(Cn/(Tmean+273))*Wspd*(es-ea))/(svp_slope+psy*(1+Cd*Wspd))
  
  # Et0 cannot be negative
  Et0 <- ifelse(Et0<0,0,Et0)
  
  
  return(Et0)
  
}
