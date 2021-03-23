import numpy as np

def rho(fw, fa, fr=None):
    """Return electrical resistivity based on fraction of water `fw`."""
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr

    r = a * rhow * phi**(-m) * (fw / phi)**(-n)
    return r

def rho_deriv_fw(fw, fa, fr):
    return rho(fw, fa, fr) * -n / fw

def rho_deriv_fr(fw, fa, fr):
    return rho(fw, fa, fr) * (n - m) / (fr - 1)

# sympy
def rho_deriv_fr_sp(fw, fa, fr):
    return rho(fw, fa, fr) * (m - n) / (1 - fr)

def sigma(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
        
    s = a**-1 * sigmaw * phi**m * (fw / phi)**n
    return s
    
def sigma0_alt(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
    
    s = (fw**n * phi**m * phi**(n-1) * sigmaw + \
         fw**(n-1) * phi**(m-1) * phi**n * rhog * (B - l) * cec) / \
         (phi**n * phi**(n-1))
    
    return s
    
def sigma0_std(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
        
    s = (fw / phi)**n * phi**m * sigmaw + \
        (fw / phi)**(n-1) * phi**(m-1) * rhog * (B - l) * cec
        
    return s
    
def sigmainf_alt(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
    
    s = (fw**n * phi**m * phi**(n-1) * sigmaw + \
         fw**(n-1) * phi**(m-1) * phi**n * rhog * B * cec) * \
         1 / (phi**(2*n-1))
    
    return s
    
def sigmainf_alt2(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
    
    s = (fw**n * phi**(m + n-1) * sigmaw + \
         fw**(n-1) * phi**(n + m-1) * rhog * B * cec) * \
         1 / (phi**(2*n-1))
    
    return s
    
def sigmainf_alt3(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
    
    s = phi**(m-n)*(fw**n * (1./rhow) + \
        fw**(n-1) * rhog * B * cec)
    
    return s
    
def sigmainf_std(fw, fa, fr=None):
    if fr is None:
        phi = fw + fa
    else:
        phi = 1 - fr
        
    s = (fw / phi)**n * phi**m * sigmaw + \
        (fw / phi)**(n-1) * phi**(m-1) * rhog * B * cec
        
    return s
    
def sigma_deriv_fw_sp(fw, fa, fr):
    return sigma(fw, fa, fr) * n / fw
    
def sigma_deriv_fr_sp(fw, fa, fr):
    return sigma(fw, fa, fr) * (n - m) / (1 - fr)

def rho_deriv_fw_alt(fw, fa, fr):
    return 1. / sigma(fw, fa, fr) * -n / fw

if __name__ == '__main__':
    rhow = 10
    sigmaw = 1./rhow
    m = 2
    n = 1.5
    a = 2
    fw = .2
    fa = .3
    fr = 1 - (fw + fa)
    
    B = 100
    l = 10
    rhog = 2
    cec = 1e3
    
    
    s0_std = sigma0_std(fw, fa, fr)
    s0_alt = sigma0_alt(fw, fa, fr)
    sinf_std = sigmainf_std(fw, fa, fr)
    sinf_alt = sigmainf_alt(fw, fa, fr)
    sinf_alt2 = sigmainf_alt2(fw, fa, fr)
    sinf_alt3 = sigmainf_alt3(fw, fa, fr)
    # ~ res = rho(fw, fa)
    # ~ res_fw = rho_deriv_fw(fw, fa, fr)
    # ~ res_fr = rho_deriv_fr_sp(fw, fa, fr)
    # ~ res_fw_alt = rho_deriv_fw_alt(fw, fa, fr)
    
    # ~ con = sigma(fw, fa)
    # ~ con_fw = sigma_deriv_fw_sp(fw, fa, fr)
    # ~ con_fr = sigma_deriv_fr_sp(fw, fa, fr)
