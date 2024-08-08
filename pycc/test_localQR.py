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

psi4.set_memory('2 GiB')
psi4.core.set_output_file('output.dat', False)
psi4.set_options({'basis': 'cc-pvdz',
                  'scf_type': 'pk',
                  'e_convergence': 1e-08,
                  'd_convergence': 1e-08,
                  'r_convergence': 1e-08
})
mol = psi4.geometry(moldict["(H2)_2"])
rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

e_conv = 1e-08
r_conv = 1e-08

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

lresp.pert_lquadresp(omega1, omega2)
lresp.lhyperpolar()
