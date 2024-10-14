"""
Test CCSD linear response functions.
"""
import numpy as np
# Import package, test suite, and other packages as needed
import psi4
from opt_einsum import contract

import pycc
#from ccwfn import ccwfn
#from lccwfn import lccwfn
#from cchbar import cchbar
#from cclambda import cclambda
#from ccdensity import ccdensity
#from ccresponse import ccresponse

import sys
#sys.path.append("/Users/jattakumi/pycc/pycc/")
sys.path.append("/Users/josemarcmadriaga/pycc_present/pycc/pycc")
from data.molecules import *
psi4.core.clean
psi4.set_memory('2 GiB')
psi4.core.set_output_file('output.dat', False)
psi4.set_options({'basis': 'cc-pvdz',
                  'scf_type': 'pk',
                  'freeze_core':'true',
                  'e_convergence': 1e-10,
                  'd_convergence': 1e-10,
                  'r_convergence': 1e-10
})
mol = psi4.geometry(moldict["(H2)_2"])
rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

e_conv = 1e-10
r_conv = 1e-10

#local
lcc = pycc.ccwfn(rhf_wfn,  local = 'PNO++', local_mos = 'BOYS', local_cutoff = 1e-07, filter=False)
lecc = lcc.lccwfn.solve_lcc(e_conv, r_conv)
lhbar = pycc.cchbar(lcc)
lcclambda = pycc.cclambda(lcc, lhbar)
llecc = lcclambda.solve_llambda(e_conv, r_conv)
ldensity = pycc.ccdensity(lcc, lcclambda)

lresp = pycc.ccresponse(ldensity)

omega1 = 0.0656
omega2 = 0.0656

lresp.pert_lquadresp(omega1, omega2, e_conv = 1e-08, r_conv = 1e-08)
lresp.lhyperpolar()

print('Time table for intermediates')
print("pertbar = %6.6f" % lresp.pertbar_t)
print("lX1 = %6.6f" % lresp.lX1_t)
print("lX2 = %6.6f" % lresp.lX2_t)
print("lY1 = %6.6f" % lresp.lY1_t)
print("lY2 = %6.6f" % lresp.lY2_t)
print("psuedoresponse = %6.6f" % lresp.pseudoresponse_t)
print("LAX = %6.6f" % lresp.LAX_t)
print("Fz_t = %6.6f" % lresp.Fz_t)
print("Bcon_t = %6.6f" % lresp.Bcon_t)
print("G_t = %6.6f" % lresp.G_t)
