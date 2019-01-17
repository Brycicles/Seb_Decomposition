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


def cx(N,control, target, control_value =1):
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

def Oracle_Decomp5qubits(choice=[], target = 4):
    '''
Using Sebastions 5 qubit decomposition. An Oracle can be setup with other
chosen values. I believe this is a phase oracle rather than a boolean oracle. 
Solution choice is made by selection of x

Needs :
    from qutip import *
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import cmath
    
    pi = math.pi 
at the top of the file
    '''
    pi = math.pi
    app = ""
    N=5
    a = 0+1j*(2*pi/8)
    b = 0+1j*(pi/8)
    dummy = np.array([[cmath.exp(a), 0],[0, cmath.exp(b)]])
    out = Qobj(dummy)
    
    for ind in range (N):
        if choice[ind] == 0:
            app = app + "x"
        else:
            app = app + "i"
        
# Components
    Choice = Gate(N,app)
    Diff1 = controlled_gate(out,N,0,4) 
    Diff2 = cnot(N,0,1) 
    Diff3 = (controlled_gate(out,N,1, 4)).conj().trans()
    Diff4 = cnot(N,0,1)
    Diff5 = cnot(N,0,2)
    Diff6 = cnot(N,1,2)
    Diff7 = controlled_gate(out,N,2,4)
    Diff8 = cnot(N,0,2)
    Diff9 = cnot(N,1,2)
    Diff10 = cnot(N,0,3)
    Diff11 = cnot(N,1,3)
    Diff12 = cnot(N,2,3)
    Diff13 = (controlled_gate(out,N,3,4)).conj().trans()
    Diff14 = cnot(N,0,3)
    Diff15 = cnot(N,1,3)
    Diff16 = cnot(N,2,3)
    Diff17 = cnot(N,1,3)
    Diff18 = cnot(N,2,3)
    Diff19 = controlled_gate(out,N,3,4)
    Diff20 = cnot(N,1,3)
    Diff21 = cnot(N,2,3)
    Diff22 = cnot(N,1,2)
    Diff23 = (controlled_gate(out,N,2,4)).conj().trans()
    Diff24 = cnot(N,1,2)
    Diff25 = controlled_gate(out,N,1,4)
    Diff26 = cnot(N,1,3)
    Diff27 = (controlled_gate(out,N,3,4)).conj().trans()
    Diff28 = cnot(N,1,3)
    Diff29 = controlled_gate(out,N,3,4)
    Diff30 = cnot(N,3,2)
    Diff31 = (controlled_gate(out,N,2,4)).conj().trans()
    Diff32 = cnot(N,3,2)
    Diff33 = controlled_gate(out,N,2,4)
    Diff34 = cnot(N,2,0)
    Diff35 = (controlled_gate(out,N,0,4)).conj().trans()
    Diff36 = cnot(N,3,0)
    Diff37 = controlled_gate(out,N,0,4)
    Diff38 = cnot(N,2,0)
    Diff39 = (controlled_gate(out,N,0,4)).conj().trans()
    Diff40 = cnot(N,1,0)
    Diff41 = controlled_gate(out,N,0,4)
    Diff42 = cnot(N,1,0)
    Diff43 = cnot(N,3,0)
    
    DiffA = Diff43*Diff42*Diff41*Diff40*Diff39*Diff38*Diff37*Diff36*Diff35
    DiffB = Diff34*Diff33*Diff32*Diff31*Diff30*Diff29*Diff28*Diff27*Diff26
    DiffC = Diff25*Diff24*Diff23*Diff22*Diff21*Diff20*Diff19*Diff18*Diff17
    DiffD = Diff16*Diff15*Diff14*Diff13*Diff12*Diff11*Diff10*Diff9*Diff8*Diff7
    DiffE = Diff6*Diff5*Diff4*Diff3*Diff2*Diff1
    
    Diff = Choice*DiffA*DiffB*DiffC*DiffD*DiffE*Choice

    return Diff

def Diffusion_Decomp5qubits(target=4, control=[0,1,2,3], control_values=[0,0,0,0]):
    '''
Using Sebastions 5 qubit decomposition.
Control[i] --> control_values[i]
Needs :
    from qutip import *
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import cmath
    
    pi = math.pi 
at the top of the file
    '''
    pi = math.pi
    flag = 0
    N=5
    app = ""
    ad = ""
    a = 0+1j*(2*pi/8)
    b = 0+1j*(pi/8)
    dummy = np.array([[cmath.exp(a), 0],[0, cmath.exp(b)]])
    out = Qobj(dummy)

    # Writing the necessary matrix to apply to accoun for control qubits.
    for chec in control:
        if control_values[chec] == 0:
            app = app + "x"
        else:
            app = app + "i"
            
    if target < all(control):
        app = "x"+app
        
    elif target > all(control):
        app = app + "x"

    else:
        print("\n Error, function only works if the target is the top qubit")
        print("\n (target=4) or the bottom qubit (target = 0)")
                
# Components
    DiffH = hadamard_transform(5)
    DiffX = Gate(N,app)
    Diff1 = controlled_gate(out,N,0,4) 
    Diff2 = cnot(N,0,1) 
    Diff3 = (controlled_gate(out,N,1, 4)).conj().trans()
    Diff4 = cnot(N,0,1)
    Diff5 = cnot(N,0,2)
    Diff6 = cnot(N,1,2)
    Diff7 = controlled_gate(out,N,2,4)
    Diff8 = cnot(N,0,2)
    Diff9 = cnot(N,1,2)
    Diff10 = cnot(N,0,3)
    Diff11 = cnot(N,1,3)
    Diff12 = cnot(N,2,3)
    Diff13 = (controlled_gate(out,N,3,4)).conj().trans()
    Diff14 = cnot(N,0,3)
    Diff15 = cnot(N,1,3)
    Diff16 = cnot(N,2,3)
    Diff17 = cnot(N,1,3)
    Diff18 = cnot(N,2,3)
    Diff19 = controlled_gate(out,N,3,4)
    Diff20 = cnot(N,1,3)
    Diff21 = cnot(N,2,3)
    Diff22 = cnot(N,1,2)
    Diff23 = (controlled_gate(out,N,2,4)).conj().trans()
    Diff24 = cnot(N,1,2)
    Diff25 = controlled_gate(out,N,1,4)
    Diff26 = cnot(N,1,3)
    Diff27 = (controlled_gate(out,N,3,4)).conj().trans()
    Diff28 = cnot(N,1,3)
    Diff29 = controlled_gate(out,N,3,4)
    Diff30 = cnot(N,3,2)
    Diff31 = (controlled_gate(out,N,2,4)).conj().trans()
    Diff32 = cnot(N,3,2)
    Diff33 = controlled_gate(out,N,2,4)
    Diff34 = cnot(N,2,0)
    Diff35 = (controlled_gate(out,N,0,4)).conj().trans()
    Diff36 = cnot(N,3,0)
    Diff37 = controlled_gate(out,N,0,4)
    Diff38 = cnot(N,2,0)
    Diff39 = (controlled_gate(out,N,0,4)).conj().trans()
    Diff40 = cnot(N,1,0)
    Diff41 = controlled_gate(out,N,0,4)
    Diff42 = cnot(N,1,0)
    Diff43 = cnot(N,3,0)
    
    DiffA = Diff43*Diff42*Diff41*Diff40*Diff39*Diff38*Diff37*Diff36*Diff35
    DiffB = Diff34*Diff33*Diff32*Diff31*Diff30*Diff29*Diff28*Diff27*Diff26
    DiffC = Diff25*Diff24*Diff23*Diff22*Diff21*Diff20*Diff19*Diff18*Diff17
    DiffD = Diff16*Diff15*Diff14*Diff13*Diff12*Diff11*Diff10*Diff9*Diff8*Diff7
    DiffE = Diff6*Diff5*Diff4*Diff3*Diff2*Diff1
    
    Diff = DiffH*DiffX*DiffA*DiffB*DiffC*DiffD*DiffE*DiffX*DiffH

    return Diff

#================================================================ Start of main
N = 5   # Also has circuit for N=4
distrib = []

# Chosen solutions
choice1 = [0,1,1,0,1]
choice2 = [1,1,0,0,0]
choice3 = [1,0,0,0,0]

# Iterate builder
# =============================================================================
Orc1 = Oracle_Decomp5qubits(choice1) # can use multiple Orcs if you want
Orc2 = Oracle_Decomp5qubits(choice2)
Orc3 = Oracle_Decomp5qubits(choice3)

Diff = Diffusion_Decomp5qubits()
# =============================================================================

#Initialise
q = qubit_states(N, np.zeros(N))
q = hadamard_transform(5)*q

# ============================================================================= Pick some initial states
# flip some states to get |00000> -> (|00101>+|00111>)/2
Set1=rx(pi/2, N, 4)
Set2=controlled_gate(snot(), N, 4,1)
Set3 = Gate(N, "xixii")
q = Set3*Set2*Set1*q
plot_fock_distribution(q)
# =============================================================================

# Iterate
#=========
for seq in range(0,20):
    q=Diff*Orc3*Orc2*Orc1*q
    plot_fock_distribution(q)
#    for ite in range(2**N):
#        pass
#        #distrib[ite] =q.full() # here get data


















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
