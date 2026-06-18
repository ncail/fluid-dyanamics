# Imports
import numpy as np

# Constants
import config as cfg

# Cross sectional area
def calculate_cross_sectional_area(d):
    r = d / 2
    A = np.pi * r**2
    return A

# Hagen-Poiseuille equation
# Describes the pressure drop across a cylindrical tube of constant cross-sectional area for laminar flow.
def calculate_hagen_poiseuille_flow(dp, r, l, mu):
    """
    Parameters:
        dp: Pressure drop across the tube (Pa)
        r: Radius of the tube (m)   
        l: Length of the tube (m)
        mu: Dynamic viscosity of the fluid (Pa*s)
    Returns:
        flow_rate: Volumetric flow rate (m^3/s)
        impedance: Hydraulic impedance (Pa*s/m^3)

    Governing relation for flow is flow Q * impedance Z = pressure drop dP
    Here, impedance Z = (8 * mu * l) / (pi * r^4) for laminar flow in a cylindrical tube.
    """
    impedance = (8 * mu * l) / (np.pi * r**4)   # Pa*s/m^3
    flow_rate = dp / impedance                  # m^3/s
    return flow_rate, impedance

# Same relation as above but the input is flow so that dp is returned.
def calculate_hagen_poiseuille_dp(flow, r, l, mu):
    """
    Parameters:
        flow: Volumetric flow rate (m^3/s)
        r: Radius of the tube (m)   
        l: Length of the tube (m)
        mu: Dynamic viscosity of the fluid (Pa*s)
    Returns:
        dp: Pressure drop across the tube (Pa)
        impedance: Hydraulic impedance (Pa*s/m^3)

    Governing relation for flow is flow Q * impedance Z = pressure drop dP
    Here, impedance Z = (8 * mu * l) / (pi * r^4) for laminar flow in a cylindrical tube.
    """
    impedance = (8 * mu * l) / (np.pi * r**4)   # Pa*s/m^3
    dp = flow * impedance                  # m^3/s
    return dp, impedance

# Reynolds number (unitless)
# This predicts whether the flow is laminar or turbulent. For laminar flow, Re < 2300. Turbulent flow occurs when Re > 4000. The range in between is transitional.
def calculate_reynolds_number(rho, v, d, mu):
    """
    The reynolds number is unitless, ensure your passed quantities cancel units.

    Parameters:
        rho: Fluid density
        v: Flow velocity (speed)
        d: Characteristic length or dimension (e.g., pipe diameter)
        mu: Dynamic viscosity
        
        A typical quantity is also nu: Kinematic viscosity (nu = mu / rho)

    For convenience, returns: 
    (
        float: Reynolds number, 
        bool : flow type
    )

    Flow type:
        True if flow is laminar
        False if turbulent
        None if transitional
    """
    reynolds = (rho * v * d) / mu

    if reynolds < 2300:
        return reynolds, True
    elif reynolds > 4000:
        return reynolds, False
    else:
        return reynolds, None  # Transitional flow

# Blasius formula, giving the Darcy friction factor for turbulent flow in smooth pipes.
# Used to calculate the Darcy–Weisbach impedance
# Accurate for Reynolds number between 4000 and 100,000 (tubulent flow).
def calculate_friction_factor(reynolds):
    factor = 0.3164 * reynolds**(-0.25)
    return factor

# Darcy–Weisbach impedance: for incompressible fluid in a smooth cylindrical pipe.
def calculate_darcy_weisbach_impedance(
    friction_factor, 
    l, 
    d, 
    rho, 
    v, 
    Q
):
    """
    Parameters:
        friction_factor: Calculated using the Blasius equation for turbulent flow calculate_friction_factor()
        l: Length of the tube (m)
        d: Diameter of the tube (m)
        rho: Density of the fluid (kg/m^3)
        v: Flow velocity (m/s)
        Q: Volumetric flow rate (m^3/s)
    
    Returns:
        impedance: (Pa*s/m3)

    Note: expression is usually proportional to 1/g.
    But since Q = rho*g, we express as proportional to rho/Q so we can input a volumetric flow rate.
    Note also that a Pascal is N/m^2 = kg/(m*s^2)
    """
    impedance = friction_factor * (l / d) * (rho * v**2) / (2 * Q)    # Pa*s/m^3
    return impedance

# Calculates pressure increase rate for a given flow rate and volume, using ideal gas law.
def pressure_increase_rate(
    flow_slpm, 
    chamber_volume_m3, 
    temperature_K=cfg.TEMPERATURE_STP, 
    molar_volume_L_per_mol=cfg.MOLAR_VOLUME
):
    """
    Calculate pressure increase rate (dP/dt) in Pascals per minute in a fixed volume chamber.

    Parameters:
    - flow_slpm: flow rate at STP (litres per minute)
    - chamber_volume_m3: chamber volume in cubic meters
    - temperature_K: temperature in Kelvin (default 298 K)
    - molar_volume_L_per_mol: molar volume at STP in litres per mole (default 22.4 L/mol)

    Returns:
    - pressure increase rate in Pascals per minute

    Notes: From ideal gas law P = (n/V) * R * T, we can derive dP/dt = (R * T / V) * dn/dt, where dn/dt is the molar flow rate.
    """

    # flow rate to moles per minute
    dn_dt = flow_slpm / molar_volume_L_per_mol  # mol/min

    # dP/dt in Pascals per minute
    dP_dt = (cfg.GAS_CONSTANT * temperature_K / chamber_volume_m3) * dn_dt  # Pa/min

    return dP_dt

# Standard Sub-Critical Gas Flow Equation. Engineering formula used to size valves for gas flow.
def std_sub_critical_gas_flow(
    Cv, 
    P1_psi, 
    P2_psi, 
    T_C, 
    G=cfg.SPECIFIC_GRAVITY
):
    """
    Calculates nitrogen gas flow (std L/min) from Cv using the average pressure method.
    Assumes flow is sub-critical (not choked) and ideal gas behavior (gas not too compressed).

    Parameters:
    Cv       : Valve flow coefficient (unitless)
    P1_psi   : Inlet pressure (psia)
    P2_psi   : Outlet pressure (psia)
    T_C      : Temperature in Celsius
    G        : Specific gravity of gas (default for N2 = 0.967)

    Returns:
    Q_scfh : Flow rate in SCFH (standard cubic feet per hour)
    dP_dependence : The term dependent on pressure drop, useful for plotting flow vs pressure drop
    """
    # Convert temperature to Rankine (absolute scale for Fahrenheit)
    T_R = (T_C + 273.15) * 1.8

    # Calculate pressure drop and average pressure
    dP = P1_psi - P2_psi
    P_avg = (P1_psi + P2_psi) / 2

    if dP <= 0:
        return 0  # No flow if backpressure is too high

    # Return term dependent on dP for convenience in plotting flow vs pressure drop
    dP_dependence = 1360 * (dP / (G * T_R))**0.5 * P_avg**0.5

    # Flow in SCFH
    # Note: 1360 comes from Hours to seconds, Pounds to ounces, Gas constants, Standard atmospheric pressure constants.
    # Also, fluid velocity is proportional to the square root of the pressure drop over the average pressure.
    Q_scfh = Cv * dP_dependence

    return Q_scfh, dP_dependence



