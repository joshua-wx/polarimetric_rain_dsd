from scipy import interpolate
import numpy as np
from pytmatrix import radar, tmatrix_aux
from pytmatrix.tmatrix import  Scatterer

def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]

#Drop shape models
def pruppacher_beard(D):
    axratio = 1.03 - 0.062 * D
    return 1 / axratio

def beard_chuang(D):
    axratio = (1.0048 + 5.7e-4*D - 2.628e-2*D**2 + 3.682e-3*D**3 - 1.677e-4*D**4)
    return 1 / axratio

def thurai_2005(D):
    axratio = 0.9707 + 4.26e-2 * D -4.29e-2 * D**2 + 6.5e-3 * D**3 - 3e-4 * D**4
    return 1/axratio

def thurai_2007(D_eq):
    """
    Axis ratio of drops with respect to their diameter.
    
    Parameter:
    ==========
        D_eq: float
            Drop diameter.
    Return:
    =======
        axratio: float
            Axis ratio of drop.
    """
    if D_eq < 0.7:
        axratio = 1.0  # Spherical
    elif D_eq < 1.5:
        axratio = 1.173 - 0.5165*D_eq + 0.4698*D_eq**2 - 0.1317*D_eq**3 - 8.5e-3*D_eq**4
    else:
        axratio = 1.065 - 6.25e-2*D_eq - 3.99e-3*D_eq**2 + 7.66e-4*D_eq**3 - 4.095e-5*D_eq**4
        
    return 1/axratio

def bzv_model(D): #good model Brandes et al. 2002
    axratio = 0.9951 + 2.51e-2 * D - 3.644e-2 * D**2 + 5.303e-3 * D**3 - 2.492e-4 * D**4
    return 1.0 / axratio

def gcb_model(D):
    axratio = 1.075 - 6.5e-2 * D - 3.6e-3 * D**2 + 4e-4 * D**3
    return 1.0 / axratio

def scatter_off_2dvd_packed(d_diameters, d_densities, scatterer):
    """
    Computing the scattering properties of homogeneous nonspherical scatterers with the T-Matrix method.
    
    Parameters:
    ===========
        d_diameters: array
            Drop diameters in mm! (or else returns values won't be with proper units.)
        d_densities: array
            Drop densities.
            
    Returns:
    ========
        dbz: array
            Horizontal reflectivity.
        zdr: array 
            Differential reflectivity.
        kdp: array
            Specific differential phase (deg/km).
        atten_spec: array
            Specific attenuation (dB/km).
    """    
    # Function interpolation.
    mypds = interpolate.interp1d(d_diameters, d_densities, bounds_error=False, fill_value=0.0)    
    scatterer.psd = mypds # GammaPSD(D0=2.0, Nw=1e3, mu=4)
    
    # Obtaining reflectivity and ZDR.
    print(scatterer)
    print(np.shape(mypds))
    print(d_diameters)
    print(d_densities)
    dbz = 10*np.log10(radar.refl(scatterer))  # in dBZ
    zdr = radar.Zdr(scatterer)  # in dB
    
    # Specific attenuation and KDP.
    scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)    
    atten_spec = radar.Ai(scatterer)  # in dB/km
    kdp = radar.Kdp(scatterer)  # in deg/km
    
    return dbz, zdr, kdp, atten_spec