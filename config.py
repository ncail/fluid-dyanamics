###################################################################
# Nitrogen gas constants
###################################################################

# Dynamic viscosity of nitrogen at STP
MU = 1.7e-5     # Pa*s, N/m^2 * s = kg/(m*s^2) * s

# Universal gas constant
GAS_CONSTANT = 8.314        # J/(mol*K)

# Molar volume of an ideal gas at STP
MOLAR_VOLUME = 22.4         # L/mol

# Temperature, STP
TEMPERATURE_STP = 273.15    # K

# Density of nitrogen at STP
RHO_STP = 1.2506            # kg/m^3

# Specific gravity of nitrogen (relative to air)
SPECIFIC_GRAVITY = 0.967

# Ratio of specific heats for nitrogen
GAMMA = 1.4     

# Molar mass of nitrogen gas
MOLAR_MASS_N2 = 28.0134    # g/mol

# Gas constant for nitrogen (R = R_universal / M) in J/(kg*K)
R_N2 = GAS_CONSTANT / (MOLAR_MASS_N2 / 1000)  # Convert g/mol to kg/mol


###################################################################
# Unit conversions
###################################################################

# Inches to meters conversion factor
inch_to_meter = 0.0254

# Feet to meters conversion factor
ft_to_m = 0.3048

# Cubic meters per second to liters per minute conversion factor
m3s_to_lpm = 60000

# PSI to Pascals conversion factor
psi_to_pa = 6894.76

# SCFH to liters per minute conversion factor
scfh_to_lpm = 28.3168 / 60

