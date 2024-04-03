import numpy as np
import matplotlib.pyplot as plt

#========== STABILITY AND CONTROL PARAMETERS ==========
S = 93.5 #wing area [m^2]
S_h = 21.72 #horizontal tail area [m^2]
A_w = 8.43 #aspect ratio of wing
A_h = 4.64 #aspect ratio of horizontal tail
S_r = S_h/S #tail volume ratio 
c = 3.80 #mean aerodynamic chord [m]
c_r = 5.38 #root chord [m]
b = 28.08 #wing span [m]
lambda_w = 0.235 #taper ratio
h_f = 3.30 #fuselage height [m]
b_f = 3.30 #fuselage width [m]
l_f = 32.50 #fuselage length [m]
b_n = 1.70*1.3 #nacelle width [m]
l_nac=5.1 #nacelle length [m]
l_n = -9.86 #nacelle arm length [m] FIND THIS VALUE!!!!
l_h = 16 #horizontal tail arm [m]
l_fn = 14.1 #distance from the wing leading to nose [m]?? not sure about this
Lambda_w = 17.5 *np.pi/180 #sweep angle [rad]
Lambda_h_4 = 26 *np.pi/180#sweep angle [rad]
Lambda_h_2 = 19.61 *np.pi/180 #sweep angle [rad]
horizontal_efficiency = 0.8 #horizontal tail efficiency
oswald_efficiency = 0.8 #horizontal tail efficiency
downwash_rate = 0 #downwash per AOA, 0 for T-tail
V_ratio = 1 #vertical tail velocity ratio
V_cruise = 213.14 #cruise speed [m/s] -f120
V_approach = 63.754 #approach speed [m/s] -f120
k_n = -2.5 #nacelle factor
T_cruise = 273.15-54.3
delta_flap = 0.5 #flap deflection [rad]
C_m_ac_flap = -0.25 *1.2 #zero lift moment coefficient due to flap deflection above (Read from graph L7 S23)
C_L_0 = 0.2 #Zero AOA lift coefficient. FIND THIS VALUE!!!!
C_L_w_cruise=0.56 # reference aerodynamic fokker
rho_0 = 1.225 #air density at sea level [kg/m^3]
rho_cruise = 0.30408 #air density at cruise altitude(35000ft) [kg/m^3]
MTOW = 361646 # Max take-off weight [N]
C_m_0_airfoil = -0.069 #FIND THIS -taken from DATCOM+assuming NACA63(2)-415
C_m_0_nacelle = 0.004+ (0.046-0.2*(c+l_n)/c)*(b_n+b_n*l_nac)/b # ASSUME LAMBDA = 0, low wing config a/c, ASSUME S_f=b_n*L_nac
C_L_h = -0.35*A_h**(1/3) #Horizontal tail lift coefficient (FORMULA FROM ADSSE L8 S17)
print(C_L_h)

a_approach = np.sqrt(1.4*287*288.15) #speed of sound at 288.15K [m/s]
a_cruise = np.sqrt(1.4*287*T_cruise) #speed of sound at 243.15K [m/s]

M_approach = V_approach/a_approach #Mach number at approach
M_cruise = V_cruise/a_cruise #Mach number at cruise

beta_cruise = np.sqrt(1-M_cruise**2) #beta at cruise
beta_approach = np.sqrt(1-M_approach**2) #beta at approach
S_overlap = b_f/np.cos(Lambda_w) *c_r #overlap area
S_net = S-S_overlap #net wing aread

ac_w_cruise = 0.40 #wing ac location at cruise (FROM GRAPH ADSEE SLIDES L7 S36)
ac_w_approach = 0.28 #wing ac location at approach

print("mach cruise:", M_cruise)
print("beta*A cruise:", beta_cruise*A_w)
print("beta*A approach:", beta_approach*A_w)
print("tan sweep/beta cruise:", np.tan(Lambda_w)/beta_cruise)
print("tan sweep/beta cruise:", np.tan(Lambda_w)/beta_approach)
print("taper wing:", lambda_w)

# ======== LOADING DIAGRAM SHIT ========

#CG VALUES
OEW = 272502 #operating empty weight [N]
OEW_kg = OEW/9.81 #operative emty weight [kg]
OEW_cg = 19.012 #operating empty weight cg location [m]
X_LEMAC = 15.5 #[m]

front_seat_cg = 0.2*l_f #x_cg of front row seats
seat_pitch = 32*0.0254 #seat pitch in meters
seat_rows = 18 #number of seat rows
passenger_weight = 90 #average passenger weight in [kg] incl. carry-on luggage

total_cargo_volume = 24.87 #total cargo volume in m^3
cargo_density = 150 #cargo density in kg/m^3

cargo_front_fraction = 0.7 #fraction of cargo in front cargo hold FIND THESE 
cargo_rear_fraction = 0.3 #fraction of cargo in rear cargo hold

cargo_weight_front = cargo_front_fraction * total_cargo_volume * cargo_density
cargo_weight_rear = cargo_rear_fraction * total_cargo_volume * cargo_density

cargo_front_cg = 0.3*l_f #cargo cg location of front cargo hold FIND THESE
cargo_rear_cg = 0.75*l_f #cargo cg location of rear cargo hold FIND THESE

fuel_mass = 4905/9.81 #Fuel mass [kg]
fuel_fraction_wing = 0 # FIND THESE FRACTIONS
fuel_fraction_center = 1

battery_mass=400 #battery mass [kg]

fuel_cg_wing = 0.55*l_f #fuel cg location in wing FIND THESE
fuel_cg_center = 23.56 #fuel cg location in center (=location of propulsion group)

battery_mass_cg = 4.6605 #removable battery cg

full_cg_list = np.array([])

def calculate_new_cg(W_original, cg_original, W_item, cg_item):
    new_cg = (W_original*cg_original+W_item*cg_item)/(W_original + W_item)
    return new_cg

def calculate_MAC_percentage(distance):
    MAC = (distance-X_LEMAC)/c
    return MAC

# ======== LOADING DUE TO CARGO ========
list_cg_cargo_front = np.array([OEW_cg ,calculate_new_cg(OEW_kg,OEW_cg,cargo_weight_front,cargo_front_cg), calculate_new_cg(OEW_kg+cargo_weight_front, calculate_new_cg(OEW_kg,OEW_cg,cargo_weight_front,cargo_front_cg), cargo_weight_rear, cargo_rear_cg)])
list_cg_cargo_rear = np.array([OEW_cg ,calculate_new_cg(OEW_kg,OEW_cg,cargo_weight_rear,cargo_rear_cg), calculate_new_cg(OEW_kg+cargo_weight_rear, calculate_new_cg(OEW_kg,OEW_cg,cargo_weight_rear,cargo_rear_cg), cargo_weight_front, cargo_front_cg)])
list_weight_cargo_front = np.array([OEW_kg,OEW_kg + cargo_weight_front, OEW_kg + cargo_weight_front + cargo_weight_rear])
list_weight_cargo_rear = np.array([OEW_kg,OEW_kg + cargo_weight_rear, OEW_kg + cargo_weight_front + cargo_weight_rear])

current_weight = OEW_kg + cargo_weight_front + cargo_weight_rear
current_cg = list_cg_cargo_front[-1]

full_cg_list = np.append(full_cg_list, list_cg_cargo_front)
full_cg_list = np.append(full_cg_list, list_cg_cargo_rear)

# ======== LOADING DUE TO PASSENGERS ========
list_cg_window_front = np.array([])
list_cg_window_back = np.array([])
list_cg_aisle_front = np.array([])
list_cg_aisle_back = np.array([])
list_cg_middle_front = np.array([])
list_cg_middle_back = np.array([])
window_weight = np.array([])
middle_weight = np.array([])
aisle_weight = np.array([])

window_weight = np.append(window_weight, current_weight)
list_cg_window_front = np.append(list_cg_window_front, current_cg)
list_cg_window_back = np.append(list_cg_window_back, current_cg)

#WINDOW SEAT LOADING
for i in range(seat_rows+1):
    cg_window_front = calculate_new_cg(current_weight, list_cg_window_front[i], passenger_weight*2, front_seat_cg + i*seat_pitch) # FRONT TO BACK multiply by 2 because 2 seats per row
    cg_window_back = calculate_new_cg(current_weight, list_cg_window_back[i], passenger_weight*2, front_seat_cg + seat_pitch*(seat_rows-i)) # BACK TO FRONT multiply by 2 because 2 seats per row
    current_weight = current_weight + passenger_weight*2 #multiply by 2 because 2 seats per row

    window_weight = np.append(window_weight, current_weight)
    list_cg_window_front = np.append(list_cg_window_front, cg_window_front)
    list_cg_window_back = np.append(list_cg_window_back, cg_window_back)

full_cg_list = np.append(full_cg_list, list_cg_window_front)
full_cg_list = np.append(full_cg_list, list_cg_window_back)

#AISLE SEAT LOADING
aisle_weight = np.append(aisle_weight, current_weight)
list_cg_aisle_front = np.append(list_cg_aisle_front, list_cg_window_back[-1])
list_cg_aisle_back = np.append(list_cg_aisle_back, list_cg_window_back[-1])
for i in range(seat_rows+1):
    cg_aisle_front = calculate_new_cg(current_weight, list_cg_aisle_front[i], passenger_weight*2, front_seat_cg + i*seat_pitch) # FRONT TO BACK multiply by 2 because 2 seats per row
    cg_aisle_back = calculate_new_cg(current_weight, list_cg_aisle_back[i], passenger_weight*2, front_seat_cg + seat_pitch*(seat_rows-i)) # BACK TO FRONT multiply by 2 because 2 seats per row
    current_weight = current_weight + passenger_weight*2 #multiply by 2 because 2 seats per row

    aisle_weight = np.append(aisle_weight, current_weight)
    list_cg_aisle_front = np.append(list_cg_aisle_front, cg_aisle_front)
    list_cg_aisle_back = np.append(list_cg_aisle_back, cg_aisle_back)

full_cg_list = np.append(full_cg_list, list_cg_aisle_front)
full_cg_list = np.append(full_cg_list, list_cg_aisle_back)

#MIDDLE SEAT LOADING
middle_weight = np.append(middle_weight, current_weight)
list_cg_middle_front = np.append(list_cg_middle_front, list_cg_aisle_back[-1])
list_cg_middle_back = np.append(list_cg_middle_back, list_cg_aisle_back[-1])


for i in range(seat_rows):
    cg_middle_front = calculate_new_cg(current_weight, list_cg_middle_front[i], passenger_weight, front_seat_cg + i*seat_pitch) # FRONT TO BACK 
    cg_middle_back = calculate_new_cg(current_weight, list_cg_middle_back[i], passenger_weight, front_seat_cg + seat_pitch*((seat_rows-1)-i)) # BACK TO FRONT 
    current_weight = current_weight + passenger_weight #only 1 seat per row. so no multiply by 2

    middle_weight = np.append(middle_weight, current_weight)
    list_cg_middle_front = np.append(list_cg_middle_front, cg_middle_front)
    list_cg_middle_back = np.append(list_cg_middle_back, cg_middle_back)

full_cg_list = np.append(full_cg_list, list_cg_middle_front)
full_cg_list = np.append(full_cg_list, list_cg_middle_back)

#======== LOADING DUE TO FUEL =========
#center tank --> wing tanks    
#fuel_mass_wing = current_weight + fuel_mass*fuel_fraction_wing
fuel_mass_center = current_weight + fuel_mass*fuel_fraction_center
fuel_mass_center_wing = current_weight + fuel_mass*(fuel_fraction_center+fuel_fraction_wing)

#cg_fuel_wing_add = calculate_new_cg(current_weight, list_cg_middle_back[-1], fuel_mass*fuel_fraction_wing, fuel_cg_wing)
cg_fuel_center_add = calculate_new_cg(current_weight,list_cg_middle_back[-1], fuel_mass*fuel_fraction_center, fuel_cg_center)

list_cg_fuel_wing_center = np.array([list_cg_middle_front[-1], cg_fuel_center_add])

fuel_weight = np.array([current_weight, fuel_mass_center_wing])

current_weight = fuel_weight[-1]
current_cg = list_cg_fuel_wing_center[-1]

#========== LOADING DUE TO BATTERY ============
battery_mass_Center = current_weight + battery_mass * battery_mass_cg

cg_battery_add = calculate_new_cg(current_weight,list_cg_fuel_wing_center[-1],battery_mass,battery_mass_cg )

list_cg_battery_center=np.array([list_cg_fuel_wing_center[-1],cg_battery_add])
battery_weight = np.array([current_weight,battery_mass_Center])

current_weight = battery_weight[-1]
current_cg = list_cg_battery_center[-1]

full_cg_list = np.append(full_cg_list, list_cg_fuel_wing_center)
full_cg_list = np.append(full_cg_list, list_cg_battery_center)
cg_min = np.min(full_cg_list)
cg_max = np.max(full_cg_list)

print("cg_min [m]:", cg_min)
print("cg_max [m]:", cg_max)

print("OEW_cg:", OEW_cg)
print("OEW_cg MAC:", calculate_MAC_percentage(OEW_cg))

print("front cargo cg:", cargo_front_cg)
print("front cargo cg MAC:", calculate_MAC_percentage(cargo_front_cg))

print("rear cargo cg:", cargo_rear_cg)
print("rear cargo cg MAC:", calculate_MAC_percentage(cargo_rear_cg))

print("fuel cg wing:", fuel_cg_wing)
print("fuel cg wing MAC:", calculate_MAC_percentage(fuel_cg_wing))
print("fuel cg center:", fuel_cg_center)
print("fuel cg center MAC:", calculate_MAC_percentage(fuel_cg_center))

#PLOTTING OF LOADING DIAGRAM
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# First plot
axs[0].plot(list_cg_cargo_front, list_weight_cargo_front, label='Cargo loading front to back')
axs[0].plot(list_cg_cargo_rear, list_weight_cargo_rear, label='Cargo loading back to front')
axs[0].plot(list_cg_window_front, window_weight, label='Window seat front to back')
axs[0].plot(list_cg_window_back, window_weight, label='Window seat back to front')
axs[0].plot(list_cg_aisle_front, aisle_weight, label='Aisle seat front to back')
axs[0].plot(list_cg_aisle_back, aisle_weight, label='Aisle seat back to front')
axs[0].plot(list_cg_middle_front, middle_weight, label='middle seat front to back')
axs[0].plot(list_cg_middle_back, middle_weight, label='middle seat back to front')
axs[0].plot(list_cg_fuel_wing_center, fuel_weight, label='Hydrogen Fuel')
axs[0].plot(list_cg_battery_center,battery_weight,label='Battery')
axs[0].axvline(x=cg_min, color='r', linestyle='--', label='CG Min')
axs[0].axvline(x=cg_max, color='b', linestyle='--', label='CG Max')
axs[0].set_xlabel('Center of gravity location [m]')
axs[0].set_ylabel('Weight [kg]')
axs[0].set_title('Loading Diagram')
axs[0].legend()

# Second plot
axs[1].plot(calculate_MAC_percentage(list_cg_cargo_front), list_weight_cargo_front, label='Cargo loading front to back')
axs[1].plot(calculate_MAC_percentage(list_cg_cargo_rear), list_weight_cargo_rear, label='Cargo loading back to front')
axs[1].plot(calculate_MAC_percentage(list_cg_window_front), window_weight, label='Window seat front to back')
axs[1].plot(calculate_MAC_percentage(list_cg_window_back), window_weight, label='Window seat back to front')
axs[1].plot(calculate_MAC_percentage(list_cg_aisle_front), aisle_weight, label='Aisle seat front to back')
axs[1].plot(calculate_MAC_percentage(list_cg_aisle_back), aisle_weight, label='Aisle seat back to front')
axs[1].plot(calculate_MAC_percentage(list_cg_middle_front), middle_weight, label='middle seat front to back')
axs[1].plot(calculate_MAC_percentage(list_cg_middle_back), middle_weight, label='middle seat back to front')
axs[1].plot(calculate_MAC_percentage(list_cg_fuel_wing_center), fuel_weight, label='Fuel')
axs[1].plot(calculate_MAC_percentage(list_cg_battery_center),battery_weight,label='Battery')
axs[1].axvline(x=calculate_MAC_percentage(cg_min), color='r', linestyle='--', label='CG Min')
axs[1].axvline(x=calculate_MAC_percentage(cg_max), color='b', linestyle='--', label='CG Max')
axs[1].set_xlabel('Center of gravity location [MAC]')
axs[1].set_ylabel('Weight [kg]')
axs[1].set_title('Loading Diagram MAC percentage')
axs[1].legend()

plt.tight_layout()
plt.show()

cg_min = calculate_MAC_percentage(cg_min) 
cg_max = calculate_MAC_percentage(cg_max)
cg = (cg_min+cg_max)/2
print("cg_min:", cg_min)
print("cg_max:", cg_max)

# ======== STABILITY FUNCTIONS ========
def calculate_C_L_alpha(A, beta, n, Lambda): #for clalpha input either wing data or tail data to get the respective clalpha
    C_L_alpha = (2*np.pi*A)/(2+np.sqrt(4+((A*beta)/n)**2*(1+(np.tan(Lambda)**2)/(beta**2))))
    return C_L_alpha

def calculate_C_L_alpha_A_h(C_L_alpha_w, b_f, b, S, S_net):
    C_L_alpha_A_h = C_L_alpha_w*(1+2.15*b_f/b)*(S_net/S)+np.pi/2*(b_f**2)/S
    return C_L_alpha_A_h

def calculate_nacelle_x_ac_chord(k_n, b_n, l_n, c, S, C_L_alpha_A_h):
    x_ac_n = 2*(k_n*(b_n**2*l_n)/(S*c*C_L_alpha_A_h)) #Multiply by 2 because of 2 times nacelle (one on each side of the fuselage)
    return x_ac_n #return X_ac/c contribution due to nacelles

def calculate_wf_x_ac_chord(C_L_alpha_A_h, sweep_angle, taper, x_ac_w, cg):
    x_ac_wf = x_ac_w - 1.8/C_L_alpha_A_h * (b_f*h_f*l_fn)/(S*c) + 0.273/(1+taper) * (b_f*cg*(b-b_f))/(c**2 * (b+2.15*b_f)) * np.tan(sweep_angle)
    return x_ac_wf #return X_ac/c contribution due to wing and fuselage

def stability_curve(C_L_alpha_h, C_L_alpha_A_h, downwash, l_h, c, tail_velocity_ratio, x_cg, x_ac, safety_factor): #S/S_h for stability, witch x = Cg/MAC
    Sh_S = 1/((C_L_alpha_h/C_L_alpha_A_h) * (1-downwash)* l_h/c * tail_velocity_ratio**2 ) * x_cg - (x_ac-safety_factor)/((C_L_alpha_h/C_L_alpha_A_h) * (1-downwash)* l_h/c * tail_velocity_ratio**2 )
    return Sh_S

# ======== STABILITY CALCULATIONS ========
C_L_alpha_w_cruise = calculate_C_L_alpha(A_w, beta_cruise, oswald_efficiency, Lambda_w)
C_L_alpha_w_approach = calculate_C_L_alpha(A_w, beta_approach, oswald_efficiency, Lambda_w)
C_L_alpha_h_cruise = calculate_C_L_alpha(A_h, beta_cruise, horizontal_efficiency, Lambda_h_2)
C_L_alpha_h_approach = calculate_C_L_alpha(A_h, beta_approach, horizontal_efficiency, Lambda_h_2)

C_L_alpha_A_h_cruise = calculate_C_L_alpha_A_h(C_L_alpha_w_cruise, b_f, b, S, S_net)
C_L_alpha_A_h_approach = calculate_C_L_alpha_A_h(C_L_alpha_w_approach, b_f, b, S, S_net)

x_ac_nacelle_cruise = calculate_nacelle_x_ac_chord(k_n, b_n, l_n, c, S, C_L_alpha_A_h_cruise)
x_ac_nacelle_approach = calculate_nacelle_x_ac_chord(k_n, b_n, l_n, c, S, C_L_alpha_A_h_approach)

x_ac_wf_cruise = calculate_wf_x_ac_chord(C_L_alpha_A_h_cruise, Lambda_w, lambda_w, ac_w_cruise, cg)
x_ac_wf_approach = calculate_wf_x_ac_chord(C_L_alpha_A_h_cruise, Lambda_w, lambda_w, ac_w_approach, cg)

x_ac_cruise = x_ac_nacelle_cruise + x_ac_wf_cruise #Total shift of ac is the shift due to nacelle + shift due to wing/fuselage
x_ac_approach = x_ac_nacelle_approach + x_ac_wf_approach

Sh_S = stability_curve(C_L_alpha_h_cruise, C_L_alpha_A_h_cruise, downwash_rate, l_h, c, V_ratio, x_ac_cruise, x_ac_cruise,0.05)

# ======== CONTROLLABILITY FUNCTIONS ========

def calculate_C_m_ac_w(C_m_0_airfoil, A, Lambda_w): #C_m_ac contribution of wing (I dont have the C_m_0_airfoil value)
    C_m_ac = C_m_0_airfoil * (A*np.cos(Lambda_w)**2)/(A+2*np.cos(Lambda_w))
    return C_m_ac

def calculate_delta_fus_C_m_ac(C_L_0, C_L_alpha_A_h, c, S): #C_m_ac contribution of fuselage
    delta_fus_C_m_ac = -1.8 * (1 - (2.5*b_f)/l_f) * (np.pi * b_f * h_f * l_f)/(4 * S * c) * C_L_0/C_L_alpha_A_h
    return delta_fus_C_m_ac 

def calculate_C_m_ac(C_m_ac_w, delta_f_C_m_ac, delta_fus_C_m_ac, delta_nac_C_m_ac):
    C_m_ac = C_m_ac_w + delta_f_C_m_ac + delta_fus_C_m_ac + delta_nac_C_m_ac
    return C_m_ac

def calculate_C_L_A_h(rho, V, S, S_h, Weight): #IDK how to find this value so i assume that C_L_(A-h) is the lift of entire aircraft minus lift of the horizontal tail and then divided by 1/2*rho*v^2*S
    L_h = C_L_h*0.5*rho*V**2*S_h #lift produced by horizontal stabilizer
    L_A_h = Weight-L_h#Lift produced by rest of aircraft
    C_L_A_h = L_A_h/(0.5*rho*V**2*S)

    return C_L_A_h

def calculate_C_m_nac(C_m_0_nac,C_l,delta_nac):
    C_m_nac=C_m_0_nac+C_l*delta_nac
    return C_m_nac


def contrallability_curve(C_L_h, C_L_A_h, lh, c, V_ratio, x_cg, x_ac, C_m_ac):
    Sh_S_contrallability = 1/(C_L_h/C_L_A_h * lh/c * V_ratio**2) * x_cg + (C_m_ac/C_L_A_h - x_ac)/(C_L_h/C_L_A_h * lh/c * V_ratio**2)
    return Sh_S_contrallability
 
# ======== CONTROLLABILITY CALCULATIONS ========

C_L_A_h_cruise = calculate_C_L_A_h(rho_cruise, V_cruise, S, S_h, MTOW)
C_L_A_h_approach = calculate_C_L_A_h(rho_0, V_approach, S, S_h, MTOW)

C_m_ac_w = calculate_C_m_ac_w(C_m_0_airfoil, A_w, Lambda_w)

C_m_ac_fuselage_cruise = calculate_delta_fus_C_m_ac(C_L_0, C_L_alpha_A_h_cruise, c, S)
C_m_ac_fuselage_approach = calculate_delta_fus_C_m_ac(C_L_0, C_L_alpha_A_h_approach, c, S)

C_m_ac_nacelle_cruise = calculate_C_m_nac(C_m_0_nacelle,C_L_w_cruise,x_ac_nacelle_cruise) #assume fusellage lift coeff is very small thus wing lift = C_L_A_h
print("Cmnac:", C_m_ac_nacelle_cruise)
C_m_ac_cruise = calculate_C_m_ac(C_m_ac_w, C_m_ac_flap, C_m_ac_fuselage_cruise, C_m_ac_nacelle_cruise)

# ======== PRINTING OF RESUTLS ========

print("C_L_alpha_w_cruise:", C_L_alpha_w_cruise)
print("C_L_alpha_w_approach:", C_L_alpha_w_approach)
print("C_L_alpha_h_cruise:", C_L_alpha_h_cruise)
print("C_L_alpha_h_approach:", C_L_alpha_h_approach)

print("C_L_alpha_A-h_cruise:", C_L_alpha_A_h_cruise)
print("C_L_alpha_A-h_approach:", C_L_alpha_A_h_approach)

print("x_ac_nacelle_cruise:", x_ac_nacelle_cruise)
print("x_ac_nacelle_approach:", x_ac_nacelle_approach)

print("x_ac_wf_cruise:", x_ac_wf_cruise)
print("x_ac_wf_approach:", x_ac_wf_approach)

print("x_ac_cruise:", x_ac_cruise)
print("x_ac_wf_approach:", x_ac_approach)


print("C_L_A-h, cruise:", C_L_A_h_cruise)
print("C_L_h", C_L_h)
print("C_l_w,cruise:", C_L_w_cruise)
print("C_l_fus,cruise:",C_L_A_h_cruise-C_L_w_cruise)
print("C_m_ac_w", C_m_ac_w)
print("C_m_ac_cruise", C_m_ac_cruise)

# ======== PLOTTING OF RESULTS ========

x_cg_values = np.linspace(-0, 1.3, 100) #Creates a range of cg values to plot the stability and controllability curves
S_Sh_values_stability_cruise = stability_curve(C_L_alpha_h_cruise, C_L_alpha_A_h_cruise, downwash_rate, l_h, c, V_ratio, x_cg_values, x_ac_cruise, 0.05)
S_Sh_values_neutral_stability_cruise = stability_curve(C_L_alpha_h_cruise, C_L_alpha_A_h_cruise, downwash_rate, l_h, c, V_ratio, x_cg_values, x_ac_cruise, 0)
S_Sh_contrallability_cruise = contrallability_curve(C_L_h, C_L_A_h_cruise, l_h, c, V_ratio, x_cg_values, x_ac_cruise, C_m_ac_cruise)

plt.plot(x_cg_values, S_Sh_contrallability_cruise, label='Controllability Cruise')
plt.plot(x_cg_values, S_Sh_values_stability_cruise, label='Stability Cruise', color='green')
plt.plot(x_cg_values, S_Sh_values_neutral_stability_cruise, label='Neutral Stability Cruise', color='yellow')
plt.fill_between(x_cg_values, S_Sh_values_stability_cruise, 0, where=(S_Sh_values_stability_cruise > 0), color='red', alpha=0.2) #Stability curve for cruise
plt.fill_between(x_cg_values, S_Sh_contrallability_cruise, 0, where=(S_Sh_contrallability_cruise > 0), color='red', alpha=0.2) #Stability curve for controllability
plt.axhline(y=S_h/S, color='grey', linestyle='--', label ='S_h/S for F100')
plt.plot([cg_min, cg_max], [S_h/S, S_h/S], color='red', label='Operative Cg range')
plt.xlabel('Center of Gravity Location [MAC]')
plt.ylabel('S_h/S')
plt.title('Stability and controllability Curve')
plt.grid(True)
plt.ylim(bottom=0)  # Limit y-axis to y > 0
plt.legend()
plt.show()

