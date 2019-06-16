"""
Spyder Editor

This is a temporary script file.
"""
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI 
import numpy as np
import scipy as sci
from scipy.interpolate import interp1d
import math 
import matplotlib  as mp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pylab
from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
from pylab import *
from scipy.optimize import curve_fit

material='Water'
T_w_in = 20 + 273
T_w_out = 33 + 273


"Part 1 - Type of Heat Esxchanger"
answer =("Shell Tube Heat Exchanger - In industry the most commonly used condenser is a water-cooled" 
         "shell and tube installed on the exhaust stream from a steam turbine in thermal" 
         "power stations. These condensers convert steam from its gaseous to its liquid state " 
         "at a pressure below atmoshperic pressure.")
print "PART 1"
print answer
print" "

"Part 2 - Determining Q_cold"
W_net = 900e6
eta = .4
#Using thermal efficiency equation
Q_hot = W_net/eta
Q_cold = (Q_hot - eta * Q_hot)
print "PART 2"
print "The rate of heat transfer in the condenser is",Q_cold*0.000001,"MW"
print " "

"Part 3 - Determine rate of heat condensation"
T_s_in = CP.PropsSI('T', 'P', 8e3, 'Q', 0.9, material)
T_s_out = CP.PropsSI('T','P',8e3,'Q',0,material)
h_g = CP.PropsSI('H','T',T_s_in,'Q',0.9,material)
h_f = CP.PropsSI('H','T',T_s_out,'Q',0,material)
m_steam = Q_cold/(h_g-h_f)
#Because the process inside the condenser is isobaric, the pressure is constant and the
#temperature remains the same during the phase change. The enthalpy of the inlet is the h_g
#and the enthalpy at the outlet is hf. To calculate the heat of condensation we can divide
#the rate of heat transfer rate by the enthalpy drop across the phase dome
print "PART 3"
print "The rate of heat condensation is", m_steam, "kg/s"
print " "

"Part 4 - Determine the mass flow rate of cooling water"
cp_water = CP.PropsSI('C','T',T_w_out,'P',101325, material)
m_water = Q_cold / (cp_water * (T_w_out-T_w_in))
T_mean_w = (T_w_in + T_w_out) / 2
#To find the mass flow rate of water, I assumed the pressure at the inlet of the tube to be 
#atmospheric pressure 
print "PART 4"
print "The mass flow rate of cooling water is", m_water, "kg/s"
print " "

"Part 5 - Determine the heat trasnsfer coefficient outside and inside the tube" 
#Assuming the dimensions of the heat exchangers, we can calculate 
#the heat trasnfer coefficients using corellations and assumptions. 

#First Guess - Based on different industry standards
D_tube_in = 0.018
D_tube_out = 0.038#1 1/2 in from the design parameters determined 
L_tube =  16.5#Values can range from 2m to 16.5m - 12 was the best value for the reduced error%
thickness = round((D_tube_out - D_tube_in)/2,4)
n_pipes = 50000
pitch = 0.048 #pitch of a 1 1/2 in outer diameter tube

#Thermal Properties of Steam and Water
#Steam:

rho_steam_in = CP.PropsSI('D', 'P', 8e3, 'Q', .9, material)
rho_steam_out = CP.PropsSI('D', 'P', 8e3, 'Q', 0, material)
viscosity_steam = CP.PropsSI('V', 'P', 8e3, 'Q', 0, material)
k_steam = CP.PropsSI('L', 'P', 8e3, 'Q', 0, material)
cp_steam = CP.PropsSI('C', 'P', 8e3, 'Q', 0, material)
Pr = cp_steam * viscosity_steam / k_steam

#Water:

rho_water = CP.PropsSI('D', 'T', T_mean_w, 'Q', 0, material)
viscosity_water = CP.PropsSI('V', 'T', T_mean_w, 'Q', 0, material)
k_water = CP.PropsSI('L', 'T', T_mean_w, 'Q', 0, material)

#Calculating Re for the inside of the tube -> Water

Re_in = 4 * m_water / (math.pi * D_tube_in * viscosity_water * n_pipes)
Pr_in = cp_water * viscosity_water / k_water

#The Reynolds number reported indicates the flow is turbulent. The nusseult correlation used
#is shown below

func = (1.82*math.log10(Re_in)-1.64)**(-2)
Nu_in = ((func/8)*(Re_in-1000)*Pr_in)/(1+12.7*(Pr_in**(2/3)-1)*(func/8)**(0.5))
h_in = Nu_in * k_water / D_tube_in
A_in = math.pi * D_tube_in * L_tube
print "PART 5"
print "The heat transfer coefficient inside one tube is", h_in, "W/m^2/K"
T_wall_in = T_mean_w + (Q_cold/(A_in * h_in * n_pipes))
print "The calculated inner wall temperature was", T_wall_in, "K"
v_in = m_water/(n_pipes*math.pi*D_tube_in**2*rho_water/4)

#After calculating the heat transfer coefficient, we can solve for the temperature of the 
#outter wall. 
#Calculating outside wall temperature using titanium pipes and the resistance of a cylinder

material_pipe = 'copper'
k_pipe = 386 #From Engineering Toolbox
Resistance_cond = math.log(D_tube_out/D_tube_in)/(2*math.pi*k_pipe*L_tube)
T_wall_out = T_wall_in + (Q_cold*Resistance_cond/(n_pipes))
print "The calculated outside wall temperature was", T_wall_out, "K"
A_out = math.pi*D_tube_out*L_tube
h_out = Q_cold/(A_out*(T_s_out - T_w_out)*n_pipes)
print "The heat transfer coefficient of condensation is", h_out, "W/m^2/K"
h_bundle = h_out*(0.6 + .42*n_pipes**(-0.25))
print "The corrected heat transfer coefficient for a bundle is", h_bundle, "W/m^2/K"
print " "
print ("Note: The inner and outter wall temperatures are very similar because of the choice"
       "made for the pipe material - very conductive - and thickness")

"Part 6 - Determine UA"
#UA was determined using LMTD method for the shell and tube heat exchanger.\
# To calculate the heat transfer coefficients, we also need the LMTD and the UA
LMTD = ((T_s_in - T_w_out) - (T_s_out - T_w_in))/math.log(
        (T_s_in - T_w_out)/(T_s_out - T_w_in))
UA = Q_cold/(LMTD)
R_in = 1/(h_in*A_in*n_pipes)
R_out = 1/(h_bundle*A_out*n_pipes)
R_conduction = Resistance_cond/n_pipes
UA_check = (R_in + R_out + R_conduction)**(-1)
error_percent = abs(UA-UA_check)*100/(UA)
print "Error:", error_percent, "%"
print " "

"Part 7 - Bundle Parameters"
#Calculated Bundle Parameters
D_shell = int(int(math.ceil(math.sqrt(n_pipes)))*pitch)
A_shell = math.pi * D_shell * L_tube
n_passes = 1

print "PART 7"
print "Number of passes:", n_passes
print "Length of tube per pass:", L_tube , "meters"
print "Footprint dimensions:", D_shell, "meters and", L_tube, "meters"
print " "

"Part 8 - pressure drop"
print "PART 8"
print ("Pressure drop across the shell side is approximately 0, because the modeled heat "
       "exchanger uses an isobaric process for the steam cycle in the condenser.")
P_drop = (0.5)*func*(L_tube/D_tube_in)*rho_water*v_in**2
print "The pressure drop across the tubes are", P_drop, "Pa"
print " "

"Part 9 - Number of tubes"
print "PART 9"
print "The number of tubes was assumed to be", n_pipes
print " "

"Part 10 - Pipe parameters"
print "PART 10"
print "The outer diameter of the pipes was determined to be", D_tube_out, "meters"
print "The inner diameter of the pipes was optimized to be", D_tube_in, "meters"
print "The thickness of the pipes were measured to be", thickness, "meters"
print "The pitch for the outter diameter difined was found to be", pitch, "meters"
print "The material for the tubes was choosen to be", material_pipe
print " "

"Part 11 - Shell parameters"
print "PART 11"
print "The diameter of the shell in a square formation was found to be", D_shell, "meters"
print " "

"Part 12 - List of thermal parameters"
print "PART 12"
print "Answers for part 12 are already given throughout the code"
print " "

"Part 13 - Identify Heat Exchangers"
print "PART 13"
print "The closest parameters to standards was found through different sources."
print "For the outter tube diameter, Design and Rating Shell and Tube Heat Exchangers - John E Edwards"
print "For the length, a 1GW GE condenser has 16.5 meters."
print "For the number of tubes, sources defined a range from 500 to 50000."
print " "

"Part 14 - Honor Code"
print "PART 14"
print ("I affirm that I did not plagiarize, use unauthorized materials, or give or recive"
       " illegitimate help on this assignment")





