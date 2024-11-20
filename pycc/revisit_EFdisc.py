#t case for computing dipole moment with simulation code"
import psi4
import numpy as np
import pycc

from data.molecules import *
import sys
sys.path.append("/Users/josemarcmadriaga/pycc_present/pycc/pycc")
import os

hf = """
F 0.0000000 0.0000000  -0.087290008493927
H 0.0000000 0.0000000 1.645494724632280
units bohr
no_reorient
symmetry c1
"""

hccf = """
C 0.0000000 0.0000000 -0.0937880
C 0.0000000 0.0000000 -1.3151600
F 0.0000000 0.0000000 1.2047720
H 0.0000000 0.0000000 -2.3892650
symmetry c1
"""

geom = "HCCF"
method = "PNO"
basis = "cc-pvdz"
init_t2 = True

e_conv = 1e-12
r_conv = 1e-12
maxiter = 2000
#0.001
#16, 0.017, 0.018, 
pert = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 
0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016,
0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 
0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035,
0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 
0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056,
0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 
0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094,
0.095, 0.096, 0.097, 0.098, 0.099,0.1]
#pert = [0.021] #, 0.086]
thresholds = [1e-07, 0] #, 0] #, 1e-09, 1e-10] #, 1e-08, 1e-09, 1e-10] #, 1e-08, 1e-09, 1e-10]
count = 0
for t in thresholds:
    PNO_cutoff = t
    for k in range(len(pert)):
        count += 1
        print("This is the pertubation", pert[k])
        ecc = []
        ecc_1 = []
        avg_pairs = []
    
        field = pert[k]
        pt_6 = [-3*field, -2*field, -field, field, 2*field, 3*field, 0] 
        for i in range(len(pt_6)):
            psi4.core.clean()
            psi4.set_memory('2 GB')
            psi4.core.set_output_file('output.dat', False)
            psi4.set_options({'basis': 'cc-pvdz',
              'scf_type': 'pk',
              'mp2_type': 'conv',
              'freeze_core': 'false',
              'e_convergence': 1e-12,
              'd_convergence': 1e-12,
              'r_convergence': 1e-12,
              'maxiter': 2000,
              'diis': 1,                
              'perturb_h': True,
              'perturb_with':'Dipole',
              'perturb_dipole': [0.0,0.0,pt_6[i]]
})
            mol = psi4.geometry(hccf)
            rhf_energy_1, rhf_wfn_1 = psi4.energy('SCF', return_wfn=True)
   
            print('This is the field', pt_6[i])

            ccsd_pos = pycc.ccwfn(rhf_wfn_1, model= 'CCSD', local_mos = None, local=method, local_cutoff=PNO_cutoff, filter=True) #, it2_opt= False) #, flag='QL', field_strength = pt_6[i])
            eccsd_pos = ccsd_pos.solve_cc(e_conv,r_conv,maxiter)

            ecc.append(eccsd_pos + rhf_energy_1)
            ecc_1.append(eccsd_pos)
            avg_pairs.append(np.average(ccsd_pos.Local.dim))
            
            print(ecc)
            print(ecc_1)
            print(avg_pairs)
#        hyperpol_tot_zz = (-ecc[0] + 8*ecc[1] -13*ecc[2] + 13*ecc[3] -8*ecc[4] + ecc[5])/(8*(field**3))
#        hyperpol_cc_zz = (-ecc_1[0] + 8*ecc_1[1] -13*ecc_1[2] + 13*ecc_1[3] -8*ecc_1[4] + ecc_1[5])/(8*(field**3))
        pairs_polar = str(geom)+"_pairs_relaxed_"+str(method)+"_"+str(PNO_cutoff)+"_"+str(basis)+"_pt6_MP2diis_None.txt"
#        tot_polar = str(geom)+"_hyperpolar_relaxed_"+str(method)+"_"+str(PNO_cutoff)+"_"+str(basis)+"_pt6_correct_1.txt"
        with open(pairs_polar, 'a') as f1:
            f1.write(str(field) + ' ')
            f1.write(str(np.average(avg_pairs)) + ' ')
            for pairs in avg_pairs:
                f1.write(str(pairs) + ' ')
            f1.write('\n')
#        with open(tot_polar,'a') as f1:
#            f1.write(str(field) + ' ')
#            f1.write(str(hyperpol_tot_zz) + ' ' + str(hyperpol_cc_zz) + '\n')
        cce_relaxed = str(geom)+"_ccenergies_relaxed_"+str(method)+"_"+str(PNO_cutoff)+"_"+str(basis)+"_pt6_MP2diis_None.txt"
        with open(cce_relaxed, 'a') as f1:
            f1.write(str(field) + ' ')
            for cc_energy in ecc_1:
                f1.write(str(cc_energy) + ' ')
            f1.write('\n')
        tote_relaxed = str(geom)+"_totenergies_relaxed_"+str(method)+"_"+str(PNO_cutoff)+"_"+str(basis)+"_pt6_MP2diis_None.txt"
        with open(tote_relaxed, 'a') as f1:
            f1.write(str(field) + ' ')
            for energy in ecc:
                f1.write(str(energy) + ' ')
            f1.write('\n')
