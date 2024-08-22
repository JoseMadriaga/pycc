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
psi4.set_options({'basis': 'aug-cc-pvdz',
                  'scf_type': 'pk',
                  'freeze_core':'true',
                  'e_convergence': 1e-12,
                  'd_convergence': 1e-12,
                  'r_convergence': 1e-12
})
mol = psi4.geometry(moldict["(H2)_2"])
#mol = psi4.geometry("""                                                 
#        O -1.5167088799 -0.0875022822  0.0744338901
#        H -0.5688047242  0.0676402012 -0.0936613229
#        H -1.9654552961  0.5753254158 -0.4692384530
#        symmetry c1
#        noreorient
#        nocom
#""")
rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

e_conv = 1e-12
r_conv = 1e-12

#local
lcc = pycc.ccwfn(rhf_wfn,  local = 'PNO', local_mos = 'BOYS', local_cutoff = 1e-05, filter=False)
lecc = lcc.lccwfn.solve_lcc(e_conv, r_conv)
lhbar = pycc.cchbar(lcc)
lcclambda = pycc.cclambda(lcc, lhbar)
llecc = lcclambda.solve_llambda(e_conv, r_conv)
ldensity = pycc.ccdensity(lcc, lcclambda)

lresp = pycc.ccresponse(ldensity)

omega1 = 0.0656
omega2 = 0.0656

lresp.pert_lquadresp(omega1, omega2, e_conv = 1e-12, r_conv = 1e-12)
lresp.lhyperpolar()
