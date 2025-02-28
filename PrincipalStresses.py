#################################################################################################################################
# This Python snipped is meant to be included in a Paraview project. When used on the grid points of a mesh, it provides the more compressional, intermediate and less compressional principal stresses
## Provided by Andres Felipe Rodriguez Corcho 
## University of Sydney
##############################################################################################################################

import numpy as np
import scipy
import builtins

dev_stress = inputs[0].CellData["projStressTensor"]

#########################################################
####### Calculation for full stress tensor with arbitrary pressure
P = np.full((len(dev_stress),),-5e8) #Arbitrary pressure to make all eigen values negative

Pf= np.array([np.full((3,3),pressure)*np.eye(3) for pressure in P])

#- Creating the full stress tensor
data_stress_tensor_full = dev_stress.copy() + Pf

#Eigen values and vectors of full tensor
eigen_val_vec_full = [np.linalg.eigh(S) for S in data_stress_tensor_full]
eigen_val_full = np.array([ev[0] for ev in eigen_val_vec_full])
eigen_vec_full = np.array([evec[1] for evec in eigen_val_vec_full])

#Eigen values and vectors of deviatoric tensor
eigen_val_vec_dev = [np.linalg.eigh(S) for S in dev_stress]
eigen_val_dev = np.array([ev[0] for ev in eigen_val_vec_dev])
eigen_vec_dev = np.array([evec[1] for evec in eigen_val_vec_dev])

#Getting sigmas for the deviatoric stress tensor -
S_mc_dev=np.copy(eigen_val_dev[:,0]) #more compressional
S_int_dev=np.copy(eigen_val_dev[:,1]) #intermediate
S_lc_dev=np.copy(eigen_val_dev[:,2]) #less compressional

#Getting sigmas for the full tensor
S_mc_full=np.copy(eigen_val_full[:,0]) #more compressional
S_int_full=np.copy(eigen_val_full[:,1]) #intermediate
S_lc_full=np.copy(eigen_val_full[:,2]) #less compressional

#Getting eigen vectors
S_mc_v=eigen_vec_dev[:,:,0]
S_int_v=eigen_vec_dev[:,:,1]
S_lc_v=eigen_vec_dev[:,:,2]

###Invariant stress in Megapascals
second_invariant_stress = 1e-6*inputs[0].CellData["projStressField"]

output.CellData.append(second_invariant_stress,"2nd Invariant stress tensor")

# #Scale mc and mc sigmas
# output.CellData.append(np.full((3,3),scale_fac), "scale mc-lc")
# #Scale intermediate sigma

S_mc,S_i,S_lc = (
        (S_mc_dev.reshape(len(S_lc_full),1)*1e-6)/second_invariant_stress*S_mc_v , 
        (S_int_dev.reshape(len(S_lc_full),1)*1e-6)/second_invariant_stress*S_int_v ,
        (S_lc_dev.reshape(len(S_lc_full),1)*1e-6)/second_invariant_stress*S_lc_v)

#Sigma saturated with the pressure
output.CellData.append(1e-6*S_mc_dev.reshape(len(S_lc_dev),1), "Dev More compressional")
output.CellData.append(1e-6*S_int_dev.reshape(len(S_lc_dev),1), "Dev Intermediate")
output.CellData.append(1e-6*S_lc_dev.reshape(len(S_lc_dev),1), "Dev Less compressional")

############################################################
#### HERE THE CALCULATED FIELDS - Sigma is exported without the pressure
#From smaller sigma to largest
output.CellData.append(S_mc, "More compressional")
output.CellData.append(S_i, "Intermediate")
output.CellData.append(S_lc, "Less compressional")


#Calculating the stress number
R = [(ev_i - ev_lc)/(ev_mc - ev_lc) for ev_mc,ev_i,ev_lc in zip(1e-6*S_mc_dev.reshape(len(S_lc_dev),1),
                                                               1e-6*S_int_dev.reshape(len(S_lc_dev),1),
                                                               1e-6*S_lc_dev.reshape(len(S_lc_dev),1))]
#print('stress_ratio is {}'.format(R))

output.CellData.append(np.array(R), "Stress number")


##Fault Regimes
FaultRegimes = np.full(np.shape(np.array(R)),3)
#Condition for Normal faulting 
#More compressional in Y must be larger than intermediate and less compressional
Nfault_ids = np.logical_and(np.abs(S_mc_v[:,1]) > np.abs(S_int_v[:,1]), 
                np.abs(S_mc_v[:,1]) > np.abs(S_lc_v[:,1]))

FaultRegimes[Nfault_ids] = 1

#Condition for strike-slip
StrikeFault_ids=np.logical_and(np.abs(S_int_v[:,1]) > np.abs(S_mc_v[:,1]), 
                np.abs(S_int_v[:,1]) > np.abs(S_lc_v[:,1]))

FaultRegimes[StrikeFault_ids] = 2

#Condition for reverse faults - whatever is not 1 or 2 conditions

output.CellData.append(FaultRegimes,"Fault regimes")


