classdef Constants
  %Collection of physical constants
  properties (Constant)
    nominal_g = 9.80665;%m/s^2
    mach_sea_level = 343;%m/s
		atm = 101.325e3;%Pa = N/m^2
		air_density_stp = 1.225;%kg/m^3 at sea level and 15 degrees C
    
		speed_of_light = 299792458;%m/s
		mu_0 = 4*pi*1e-7;%Henry/m
		ep_0 = 1/(Constants.mu_0*Constants.speed_of_light^2);%Farad/m
		
		Gravitional = 6.6740831e-11;%m^3/(kg*s^2) Gravitaional constant  - F=G*m1*m2/r^2
		
		q0 = 1.60217662e-19;%Coulombs
		avagadro = 6.02214076e23;%1/mol
		boltzmann_constant = 1.38064852e-23;%J/K
    standard_temperature = 290;%K
	end
end