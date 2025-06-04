#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
get_ipython().run_line_magic('matplotlib', 'inline')

from qiskit import QuantumCircuit, transpile


from qiskit_aer import AerSimulator


from qiskit.visualization import plot_histogram
from qiskit.quantum_info import Statevector
import matplotlib as mpl







def qberry(sites,D):
    H=0
    Hsolo = 0
    Sx=(1/np.sqrt(2))*np.array([[0,1,0],[1,0,1],[0,1,0]])
    Sy=(1/np.sqrt(2))*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
    Sz=np.array([[1,0,0],[0,0,0],[0,0,-1]])
    I=np.array([[1,0,0],[0,1,0],[0,0,1]])
    Sup = np.sqrt(2)*np.array([[0,1,0],[0,0,1],[0,0,0]])
    Sdown = np.sqrt(2)*np.array([[0,0,0],[1,0,0],[0,1,0]])
    #nearest neighbor interaction XYZ
    output = np.array(np.zeros(3**sites),dtype=complex)

 ##########_________________________________________Define SxSy   -----------------------------------------------
   
    for i in range(sites-1):
        Hx = np.kron(Sx,Sx)
        Hy = np.kron(Sy,Sy)
        Hz = np.kron(Sz,Sz)
        #left kron
        if i == 0:
            pass
        else:
            for j in range(i):    
                Hx = np.kron(I,Hx)
                Hy = np.kron(I,Hy)
                Hz = np.kron(I,Hz)                
        #right
        if i==(sites-2):
            pass
        else:
            for j in range(sites-2-i):
                Hx = np.kron(Hx,I)
                Hy = np.kron(Hy,I)
                Hz = np.kron(Hz,I)
               
               
        H = H+Hx+Hy+Hz

   
               
       
    #B is variable, goes from 0 to 2 in the 2 site spin half example
       
#######_______________________________________________________________________Define Sz^2  ---------------------------  
    for i in range(sites):
        Dterm = Sz
        #left
        if i==0:
            pass
        else:
            for j in range(i):
                Dterm = np.kron(I,Dterm)
       
        #right
        if i==(sites-1):
            pass
        else:
            for j in range(sites-1-i):
                Dterm = np.kron(Dterm,I)
               
        Hsolo = Hsolo+Dterm**2
   
   
   
   
####______________________________________________________Define berry phase bond-------------------------------------
   
    plusterm = np.kron(Sup,I)
    minusterm = np.kron(Sdown,I)
    zterm = np.kron(Sz,I)

    for i in range(sites-3):
        plusterm = np.kron(plusterm,I)
        minusterm = np.kron(minusterm,I)
        zterm = np.kron(zterm,I)

    plusterm = np.kron(plusterm,Sdown)
    minusterm = np.kron(minusterm,Sup)
    zterm = np.kron(zterm,Sz)

##________________________________________________Diagonalization using different hamiltonian                
       

    

    
    H81 = 0.5*((minusterm) + (plusterm))+zterm
    Ham=H+H81+D*Hsolo


    w,v = np.linalg.eig(Ham)
    

    wmin=np.min(w)
    ind=np.argmin(w)
    vmin=v[:,ind]
       
       
    result=vmin   
    return result, wmin





# In[5]:


def purity(matrix):

    #purity part
    

    Nsite = 3
    Nq = 17 

    #11 = 2*5 + 1

    qc =  QuantumCircuit(Nq,1)

    state1 = np.zeros(2**8,dtype=complex)
    state2 = np.zeros(2**8,dtype=complex)




    state1[0:243] = matrix
    state2[0:243] = matrix



    qc.initialize(state1, [i for i in range(0,8)])
    qc.initialize(state2, [i for i in range(8,16)])


    qc.h(16)

    qc.barrier()

    for i in range(8):
        qc.cswap(16,i,i+8)

    qc.barrier() 

    qc.h(16)
    
    
    #backend = Aer.get_backend('statevector_simulator')
    #results = execute(qc,backend).result()


    s = Statevector(qc)

    #s = results.get_statevector()
    a = s.probabilities([10])[0]
    b = s.probabilities([10])[1]
    return a-b


    #style = {'fontsize':20, 
     #        'subfontsize':10,
     #       }

    #qc.draw(output='mpl', fold=50, style=style)


# In[37]:


def dephasedpurity(matrix):
    



    Nsite = 5
    Nq = 18 

    #21 = 4*5 + 1

    qc =  QuantumCircuit(Nq,1)

    state1 = np.zeros(2**8,dtype=complex)
    state2 = np.zeros(2**8,dtype=complex)




    state1[0:243] = matrix
    state2[0:243] = matrix



    qc.initialize(state1, [i for i in range(0,8)])
    qc.initialize(state2, [i for i in range(8,16)])


    qc.h(17)
    qc.barrier()

    for i in range(16):

        qc.cx(i,16)
        qc.reset(16)

    qc.barrier()

    for i in range(8):
        qc.cswap(17,i,i+8)

    qc.barrier() 

    qc.h(17)
    qc.measure(17,0)

    simulator = AerSimulator()
    circ = transpile(qc, simulator)
    
    # Run and get counts
    shots=1000
    result = simulator.run(circ,shots=shots).result()
    counts = result.get_counts(circ)
    output = (counts['0']-counts['1'])/shots
    return output   

    


# In[42]:


datapoint = 80
x = np.linspace(0.5,1.1,datapoint)
coherence = np.zeros(datapoint)
energy = np.zeros(datapoint)

index = 0
for i in x:
    
    matrix, e = qberry(5,i)  #call the main function
    p = 1 - dephasedpurity(matrix)
    coherence[index] = p
    energy[index] = e
    index=index+1

    


# In[45]:


fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, figsize=(8,16))

ax1.plot(x, coherence)
ax1.set_title('L=5')
ax1.set_xlabel('D', fontsize = 20)
ax1.set_ylabel(r'$C_{l2}$', fontsize = 20)


ax2.plot(x, -np.gradient(coherence))
ax2.set_title('L=5')
ax2.set_xlabel('D', fontsize = 20)
ax2.set_ylabel(r'-$\frac{\partial C_{l2}}{\partial D}$', fontsize = 20)


ax3.plot(x, -np.gradient(np.gradient(energy)))
ax3.set_title('L=5')
ax3.set_xlabel('D', fontsize = 20)
ax3.set_ylabel(r'-$\frac{\partial^{2} E}{\partial D^{2}}$', fontsize = 20)

ax4.plot(x, -np.gradient(np.gradient(coherence)))
ax4.set_title('L=5')
ax4.set_xlabel('D', fontsize = 20)
ax4.set_ylabel(r'-$\frac{\partial^{2} C_{l2}}{\partial D^{2}}$', fontsize = 20)


plt.subplots_adjust(hspace=0.8,left=0.2,
                    bottom=0.2, 
                    right=0.9, 
                    top=0.9)
fig.savefig('2025data4.png', facecolor='w')


# In[ ]:




