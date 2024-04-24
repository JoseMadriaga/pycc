"""
Test simulation code of CCSD quadratic response function.
"""

# Import package, test suite, and other packages as needed
import psi4
import pycc
import sys 
sys.path.append("/Users/josemarcmadriaga/pycc_present/pycc/pycc")
from data.molecules import * 

  
#geom = []
#local_model = []
threshold = [1e-03, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 0]

geom = "(H2)_4"
lmo = "PNO++"

psi4.set_memory('2 GiB')
psi4.core.set_output_file('output.dat', False)
psi4.set_options({'basis': 'aug-cc-pvdz',
                  'scf_type': 'pk',
                  'freeze_core': 'true',
                  'e_convergence': 1e-08,
                  'd_convergence': 1e-08,
                  'r_convergence': 1e-08})
mol = psi4.geometry(moldict[geom])
rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

for t in threshold:
    e_conv = 1e-08
    r_conv = 1e-08
    maxit = 1000

    #OR frequencies
    omega1 = -0.0656
    omega2 = 0.0656

    cc = pycc.ccwfn(rhf_wfn, local_mos = 'BOYS', local= lmo, local_cutoff = t, filter=True, omega = omega2)
    ecc = cc.solve_cc(e_conv, r_conv, maxit)
    hbar = pycc.cchbar(cc)
    cclambda = pycc.cclambda(cc, hbar)
    lecc = cclambda.solve_lambda(e_conv, r_conv, maxit)
    density = pycc.ccdensity(cc, cclambda)
    
    resp = pycc.ccresponse(density)
    
    resp.pert_quadresp(omega1, omega2, e_conv, r_conv, maxit)
    OR = resp.hyperpolar()
   
    B_avg = str(geom)+"_omega1_adz_OR_"+str(lmo)+".txt"

    with open(B_avg, 'a') as f1:
        f1.write(str(t) + ' '+ str(cc.Local.T2_ratio)+' ') 
        
        for i in range(0,3):
             if i != 2:
                 f1.write(str(OR[0][2,i,i]) + ' ' +str(OR[0][i,2,i]) +' '+str(OR[0][i,i,2])+' ')  
        f1.write(str(OR[0][2,2,2]) + ' '+str(OR[1])+ '\n')  
     
    del cc, hbar, cclambda, density, resp  
