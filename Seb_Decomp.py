# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Creates a random NxN density matrix.
rand_dm(N, density=0.75, pure=False, dims=None)

Creates a random NxN sparse unitary quantum object
rand_unitary(N, density=0.75, dims=None)

sesolve(H, rho0, tlist, e_ops=[], args={}, options=None, progress_bar=<qutip.ui.progressbar.BaseProgressBar object>, _safe_mode=True)

Hadamard gate is snot() (Don't use the arguments for snot. Snot(N) not sure what it does. For large hadamard gayes use tensor (snot(),snot()))

gate_sequence_product(U_list, left_to_right=True)
Calculate the overall unitary matrix for a given list of unitary operations
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

np.set_printoptions(threshold=np.nan)

pi = math.pi 


def V(N, control, target, control_value=0):    
    a = 0+1j*(2*pi/8)
    b = 0+1j*(pi/8)
    dummy = np.array([[cmath.exp(a), 0],[0, cmath.exp(b)]])
    out = Qobj(dummy)
    
    q = controlled_gate(out, N,control, target, control_value)
    
    return q


def V2(N, control, target, control_value=0):    
    a = 0+1j*(2*pi/4)
    b = 0+1j*(pi/4)
    dummy = np.array([[cmath.exp(a), 0],[0, cmath.exp(b)]])
    out = Qobj(dummy)
    
    q = controlled_gate(out, N,control, target, control_value)
    
    return q

def cx(N,control, target, control_value =0):
    q = controlled_gate(sigmax(), N,control, target, control_value)   
    return q

def Gate(N, sequence = [] ):
    '''
Input: N = number of qubits in the register.
Input: sequence = string of characters to label the gates
    example: sequence = "xii" would apply x gate to qubit 1 and identity gate 
                        to idenity gates to qubits 2 and 3. 
                        |000> = sequence(|100>) 
                        Also using my cheap example, the qubit on top is the rightmost qubit
    '''
    if len(sequence) != N:
        print("\n\t\t Error! The number of gates isn't matching the sequence.")
        print("\t\t See the function 'Gate' in the script.")

    elif N<2:
        print("\n Need to have a gate spanning more than 1 qubit.")

    else:
        if sequence[0] == 'x':
            q = sigmax()
        elif sequence[0] == 'y':
            q = sigmay()
        elif sequence[0] == 'z':
            q = sigmaz()
        elif sequence[0] == 'h':
            q = snot()
        elif sequence[0] == 'i':
            q = qeye(2)

        for i in range (1,N):
            if sequence[i] == 'x':
                q = tensor(q, sigmax())

            elif sequence[i] == 'y':
                q = tensor(q, sigmay())

            elif sequence[i] == 'z':
                q = tensor(q, sigmaz())

            elif sequence[i] == 'h':
                q = tensor(q, snot())

            elif sequence[i] == 'i':
                q = tensor(q, qeye(2))

    return q


#==============
N = 5   # Also has circuit for N=4

#Initialise
q = qubit_states(N, np.zeros(N))
# flip some states to get |00000> -> (|00101>+|00111>)/2
print(q)

Set1=rx(pi/2, N, 4)
Set2=controlled_gate(snot(), N, 4,1)
Set3 = Gate(N, "xixii")
q = Set3*Set2*Set1*q
plot_fock_distribution(q)


# Diffusion operator coonstruct (4 qubits)
# ===========================
DiffH = hadamard_transform(N)

DiffX = Gate(N, "xiiii")


Diff1 = V(N,0,4) 
Diff2 = cx(N,0,1) 
Diff3 = (V(N,1, 4)).conj().trans()
Diff4 = cx(N,0,1)
Diff5 = cx(N,0,2)
Diff6 = cx(N,1,2)
Diff7 = V(N,2,4)
Diff8 = cx(N,0,2)
Diff9 = cx(N,1,2)
Diff10 = cx(N,0,3)
Diff11 = cx(N,1,3)
Diff12 = cx(N,2,3)
Diff13 = (V(N,3,4)).conj().trans()
Diff14 = cx(N,0,3)
Diff15 = cx(N,1,3)
Diff16 = cx(N,2,3)
Diff17 = cx(N,1,3)
Diff18 = cx(N,2,3)
Diff19 = V(N,3,4)
Diff20 = cx(N,1,3)
Diff21 = cx(N,2,3)
Diff22 = cx(N,1,2)
Diff23 = (V(N,2,4)).conj().trans()
Diff24 = cx(N,1,2)
Diff25 = V(N,1,4)
Diff26 = cx(N,1,3)
Diff27 = (V(N,3,4)).conj().trans()
Diff28 = cx(N,1,3)
Diff29 = V(N,3,4)
Diff30 = cx(N,3,2)
Diff31 = (V(N,2,4)).conj().trans()
Diff32 = cx(N,3,2)
Diff33 = V(N,2,4)
Diff34 = cx(N,2,0)
Diff35 = (V(N,0,4)).conj().trans()
Diff36 = cx(N,3,0)
Diff37 = V(N,0,4)
Diff38 = cx(N,2,0)
Diff39 = (V(N,0,4)).conj().trans()
Diff40 = cx(N,1,0)
Diff41 = V(N,0,4)
Diff42 = cx(N,1,0)
Diff43 = cx(N,3,0)

DiffA = Diff43*Diff42*Diff41*Diff40*Diff39*Diff38*Diff37*Diff36*Diff35*Diff34*Diff33*Diff32
DiffB = Diff31*Diff30*Diff29*Diff28*Diff27*Diff26*Diff25*Diff24*Diff23*Diff22
DiffC = Diff21*Diff20*Diff19*Diff18*Diff17*Diff16*Diff15*Diff14*Diff13*Diff12
DiffD = Diff11*Diff10*Diff9*Diff8*Diff7*Diff6*Diff5*Diff4*Diff3*Diff2*Diff1

Diff = DiffH*DiffX*DiffA*DiffB*DiffC*DiffD*DiffX*DiffH


# =========================works backwards=====================================
# #    ADiff= Diff1*Diff2*Diff3*Diff4*Diff5*Diff6*Diff7*Diff8*Diff9*Diff10*Diff11
# #    BDiff= Diff12*Diff13*Diff14*Diff15*Diff16*Diff17*Diff18*Diff19*Diff20*Diff21
# #    CDiff= Diff22*Diff23*Diff24*Diff25*Diff26*Diff27*Diff28*Diff29*Diff30*Diff31
# #    DDiff= Diff32*Diff33*Diff34*Diff35*Diff36*Diff37*Diff38*Diff39*Diff40*Diff41*Diff42*Diff43
# #    
# #    fidd = DiffH*DiffX*ADiff*BDiff*CDiff*DDiff*DiffX*DiffH
# =============================================================================
matrix_histogram_complex(Diff)

## Loop
##=========
#for seq in range(0,10):
q=Diff*q
plot_fock_distribution(q)

# =============================================================================
# # This is the 4 qubit setup
# N=4
# DiffX = Gate(N, "xxxx")
# 
# Diff1 = V2(N, 0, 3)
# Diff2 =  cnot(N, control = 0, target = 1)
# Diff3 = (V2(N, 1, 3)).conj().trans()
# Diff4 = cnot(N, control = 0, target = 1)
# Diff5 = V2(N,1,3)
# Diff6 =  cnot(N, control = 1, target = 2)
# Diff7 = (V2(N,2,3)).conj().trans()
# Diff8 =cnot(N, control = 0, target = 2)
# Diff9 = V2(N,2,3)
# Diff10 = cnot(N, control = 1, target = 2)
# Diff11 = (V2(N,2,3)).conj().trans()
# Diff12 =  cnot(N, control = 0, target = 2)
# Diff13 = V2(N,2,3)
# 
# Diffa2= Diff13*Diff12*Diff11*Diff10*Diff9*Diff8*Diff7*Diff6*Diff5*Diff4*Diff3*Diff2*Diff1
# 
# Diff1 = V2(N, 0, 3)
# Diff2 = cx(N,0,1) # cnot(N, control = 0, target = 1)
# Diff3 = (V2(N, 1, 3)).conj().trans()
# Diff4 = cx(N,0,1) # cnot(N, control = 0, target = 1)
# Diff5 = V2(N,1,3)
# Diff6 = cx(N,1,2) # cnot(N, control = 1, target = 2)
# Diff7 = (V2(N,2,3)).conj().trans()
# Diff8 = cx(N,0,2) # cnot(N, control = 0, target = 2)
# Diff9 = V2(N,2,3)
# Diff10 = cx(N,1,2) # cnot(N, control = 1, target = 2)
# Diff11 = (V2(N,2,3)).conj().trans()
# Diff12 = cx(N,0,2) # cnot(N, control = 0, target = 2)
# Diff13 = V2(N,2,3)
# 
# #
# Diffa3 = DiffX*Diff13*Diff12*Diff11*Diff10*Diff9*Diff8*Diff7*Diff6*Diff5*Diff4*Diff3*Diff2*Diff1*DiffX
# 
# =============================================================================
