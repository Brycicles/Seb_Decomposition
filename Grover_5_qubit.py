# -*- coding: utf-8 -*-
"""
Spyder Editor



"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import random
import time

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
 #   print(out)
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

def ConvertChoicetoDec(choice):
    '''
Specifically needs to be a binnary number of 5 digits, as that is relevant here
Does so by elementwise multiplation the summing the product
    '''
    BinRepr = [16,8,4,2,1]
    dummy = np.multiply(choice, BinRepr)
    Decimal = np.sum(dummy)
    
    return Decimal

def RandDistrib(Num_of_qubits, strongWeight=[], weakWeight=[]):
    '''
Produces an arbitrary distribution. Has option to give some some states a
higher or lower probability than a uniform distrib between 0 and 1.
This is done with the weights which allow an if sub condition that is a number
appears in their list then the uniform distribution will be between 0->0.3 and 
0.7 to 1.0 for the weak and strong weights respectively. Can change of course.
    '''
    Num_of_states = 2**N
    Register = []
    for state in range(0,Num_of_states):
        if random.uniform(0,1) < 0.5:  # give some negative amplitudes
            flip = -1
        else:
            flip = 1
            
        if any(strongWeight) == state:
            Register = np.append(Register,[flip*(random.uniform(0.7,1))])
        elif any(weakWeight) == state:
            Register = np.append(Register, [flip*(random.uniform(0,0.3))])
        else:
            Register = np.append(Register, [flip*(random.uniform(0.3,0.7))])
        print("\n\n At state ", state, " Register is ", Register)
            
    norm = np.dot(Register,(np.transpose(Register)))
    Register = np.multiply(Register,(1/(math.sqrt(norm))))
    probs = np.multiply(Register,Register)
    if sum(probs)!= 1.0:
        print("\n Error, RandDistr didn't output normalised state.")
        print("\n Sum of probs ", sum(probs))
    
    Register = np.transpose(Register)        
    print("\n Register", Register)
    
    return Register

def Qubit_Builder(var = 0.1, N=5):
    '''
Builds a arbitrary distribution. N = number of qubits in register (fixed here.
var is the variance allowed around uniform state. Found using a random matrix
of [  rand ,  0  ]  *  [1] =  [ rand ]
   [ 1-rand,  0  ]  *  [0] =  [1-rand]
    '''
    qu=[]
    for qub in range(0,N):
        rando = random.uniform(0.5-var,0.5+var) # set to prevent some dramatic differences
        dummy = Qobj([[np.sqrt(rando),0], [np.sqrt(1-rando),0]])*basis(2,0)
        qu = np.append(qu, dummy)
        print(rando)
        
    q = qu[0]
    for step in range(1,N):
        q = tensor(q,qu[step])
        
 #   if q.dag()*q != basis(1,0):
    if abs((q.dag()*q)[0][0][0]) != 1:
        print("\n Warning! The function Qubit_Builder() did not out put a ")
        print("\n normalised state. q.dag()*q = ", (q.dag()*q))
        
    print(q)    
    return q

def User_input(q, qubits):
    '''
Function to recieve choice of multiple solutions. Allows automation to find ave
of marked and unmarked states for any number of solns. 
q is the quantum state from initial distrib. qubits however is the number of
qubits used in the register
k_dex is list of the marked states, and k is the amplitudes for these states
l_dex is list of the unmarked states, and l is the amplitude for these states
    '''
    print("\n Give the marked solutions. Enter in the following example:")
    print("\n Choice n: \n X X X X X")
    print("\n Where n is the number of choices made and XXXXX is any")
    print("\n combination of X=0 or 1. Must use 5 X for each choice. If you")
    print("\n have given all the solutions, write 'end' in place of 'X X X X X'")
    print("\n\n Alternatively, to use the default setting in the function,")
    print("\n enter: default")
    n = 0
    default = [0,0,1,1,1,1,1,0,1,0,0,1,0,1,1,0,0,1,0,1,1,1,0,0,1] 
    Default_number = 5 # 5 solns
    
    stuck = 'n'
    err_counter = 0
    choice = []
    while True:
        if err_counter >1:
            print("\n Do you want to end? y/n")
            stuck = input()
            if stuck == 'y':
                break
            else:
                err_counter = 0
            
        print("\n Choice ", n+1, ":")
        inp = input()
        if inp == 'end':
            break
        if inp == 'default':
            choice = default
            n=Default_number
            break
        
        inp = list(map(int, inp.split(' ')))
        if len(inp)!= 5:
            print("\n Error in input. Must input 5 values, each seperated by")
            print("\n a space. For example, \nChoice x: \n 0 1 1 1 0")
            err_counter = err_counter + 1
            continue
        
        choice = np.append(choice,inp)
        n=n+1
        
    if n == 0:
        print("\n No solutions given. Error!")

    print("\n ",n," solutions recieved:")
    i=0
    k_dex = []
    k = []
    l_dex = []
    l = []
    for pr in range (0,n):
        print("\n", choice[i:i+5])
        State = ConvertChoicetoDec(choice[i:i+5])
        k_dex.append(int(State))
        k.append(abs(q[State][0][0]))
        i = i+5
        
    for seq in range(0,(2**qubits)):
        l_dex.append(seq)
    for wq in range(0,len(k_dex)):
        l_dex.remove(k_dex[wq])        
    for seq in range(0,len(l_dex)):
        l.append(abs(q[seq][0][0]))

    #collect averages
    l_ave = sum(l)/len(l)
    k_ave = sum(k)/len(k)
    
    #collect variances
    sigk = []
    sigl = []
    for al in range(0,len(k_dex)):
        sigk.append(abs(k[al]-k_ave)**2)
    for al in range(0,len(l_dex)):
        sigl.append(abs(l[al]-l_ave)**2)

    sigmak = (1/n)*sum(sigk)
    sigmal = (1/((2**qubits)-n))*sum(sigl)
            
    return choice, k_dex,l_dex, n, k_ave, l_ave, sigmak, sigmal
        
def Oracle(choice, n, Num):
    '''
    '''
    i=0
    Orc = qeye(2)
    for ident in range(0,Num-1):  # need for compatible shapes
        Orc=tensor(Orc,qeye(2))
    for seq in range(0,n):
        Orc = Oracle_Decomp5qubits(choice[i:i+5])*Orc
        i = i+5
    
    return Orc

#================================================================ Start of main
N = 5
n = 2**N
distrib = []
Nume = 20

Diff = Diffusion_Decomp5qubits()
# =============================================================================
#Initialise
q = Qubit_Builder()
plot_fock_distribution(q)

choice,mk_soln,unmk_soln,num_solns,k_ave,l_ave,sigmak,sigmal = User_input(q, N)
Data_amplitude = np.transpose((q.full()).real)
Data_amplitude_ini = np.transpose((q.full()).real)[0]

# Predicted optimal number of states and max probability
# (from the paper, eq.25 and 21 :arXiv:quant-ph/9801066v2 11 May 1998)
# ========================================
T1=pi/2#-0.5*(k_ave/l_ave)
T2=-np.arctan((k_ave/l_ave)*np.sqrt(num_solns/(n-num_solns)))#(pi/4)*(np.sqrt(n/num_solns))
T3=np.arccos(1-2*(num_solns/n))#-(pi/24)*(np.sqrt(num_solns/n))
##T4=np.sqrt(n/num_solns)
T=(T1+T2)/T3#T1+T2+T3#+T4

Pmax_estimate = 1 - ((n)-num_solns)*(sigmak**2)

print("\n The expected number of iterations needed to peak in probability is")
print("\n ", T, " and the max probability (summed from each marked soln) is")
print("\n ", Pmax_estimate)

# =============================================================================
# Oracle         
Orc = Oracle(choice, num_solns, N)

# Iterate
#==============
peak = 1
old = 0

for seq in range(0,Nume):
    q=Diff*Orc*q
    Data_amplitude = np.vstack((Data_amplitude,np.transpose((q.full()).real)))
    plot_fock_distribution(q)
    new = np.absolute(q[mk_soln[0]])
    if new < old and peak < 3:
        print("\n peak ", peak, " seen at seq", seq)
        Data_amplitude_fin = np.transpose((q.full()).real)[0]
        Data_amp_change = Data_amplitude_fin-Data_amplitude_ini
        Data_prob_change1 = np.absolute(Data_amplitude_fin) - np.absolute(Data_amplitude_ini)
        Data_prob_change2 = np.absolute(Data_amp_change)
        peak = peak + 1
        
    old = new

np.savetxt('results_amp.csv', Data_amplitude, delimiter=',')
Data_probability = np.absolute(Data_amplitude)
np.savetxt('results_prob.csv', Data_probability, delimiter=',')
np.savetxt('Data_amp_change.csv' , Data_amp_change, delimiter=',')
np.savetxt('Data_prob_change1.csv' , Data_prob_change1, delimiter=',')
np.savetxt('Data_prob_change2.csv' , Data_prob_change2, delimiter=',')

# =============================================================================
# # Plotting 
# # =============================================================================
# #plt.gcf().clear()
# 
# for coun in range(0,n-num_solns):
#     plt.plot(list(range(0,Nume+1)),Data_probability.transpose()[unmk_soln[coun]],'b--',label='Unmarked #' +str(coun))    
#     
# for coun2 in range(0,num_solns):
#     plt.plot(list(range(0,Nume+1)),Data_probability.transpose()[mk_soln[coun2]],'r--',label='Marked #'+str(coun2))
#     
# plt.axis([0, n, 0, 1])
# plt.legend(loc='upper right')
# plt.ylabel('Probability')
# plt.xlabel('Iterate)')
# plt.show()
# =============================================================================

# =============================================================================
