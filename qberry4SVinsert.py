#!/usr/bin/env python
# coding: utf-8

# In[1]:


from qiskit import QuantumCircuit, transpile
from qiskit import Aer, execute
from qiskit.visualization import plot_histogram
from qiskit.providers.aer import QasmSimulator
from qiskit.quantum_info import Statevector
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import cmath


# In[2]:





def qberry(sites,steps,D):
    H=0
    Hsolo = 0
    Sx=(1/np.sqrt(2))*np.array([[0,1,0],[1,0,1],[0,1,0]])
    Sy=(1/np.sqrt(2))*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
    Sz=np.array([[1,0,0],[0,0,0],[0,0,-1]])
    I=np.array([[1,0,0],[0,1,0],[0,0,1]])
    Sup = np.sqrt(2)*np.array([[0,1,0],[0,0,1],[0,0,0]])
    Sdown = np.sqrt(2)*np.array([[0,0,0],[1,0,0],[0,1,0]])
    #nearest neighbor interaction XYZ
    output = np.array(np.zeros((steps,3**sites)),dtype=complex)

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
       

    for i in range(steps):

        phi = np.pi*i*2/(steps-1)
        #phi introduces the berry phase, for 4 sites, we have 0, 2pi/3, 4pi/3, 2pi
        H81 = 0.5*((np.exp(1j*phi)*minusterm) + (np.exp(-1j*phi)*plusterm))+zterm
        Ham=H+H81+D*Hsolo


        w,v = np.linalg.eig(Ham)
        #spa=csr_matrix(Ham)
        #w, v = eigs(spa,k=6,which='SR')

        wmin=np.min(w)
        ind=np.argmin(w)
        vmin=v[:,ind]
        output[i,:] = vmin
       
       
    result=output    
    return result




# In[3]:


def recircuit():

    #from qiskit import QuantumCircuit, transpile, IBMQ
    #from qiskit import Aer, execute
    #from qiskit.visualization import plot_histogram
    #from qiskit.providers.aer import QasmSimulator
    #from qiskit.quantum_info import Statevector
    #import matplotlib as mpl

    Nsite = 4
    Nq = 29
    qc =  QuantumCircuit(Nq,1)

    state1 = np.zeros(2**7,dtype=complex)
    state2 = np.zeros(2**7,dtype=complex)
    state3 = np.zeros(2**7,dtype=complex)
    state4 = np.zeros(2**7,dtype=complex)



    state1[0:81] = matrix[0]
    state2[0:81] = matrix[1]
    state3[0:81] = matrix[2]
    state4[0:81] = matrix[3]


    qc.initialize(state1, [i for i in range(0,7)])
    qc.initialize(state2, [i for i in range(7,14)])
    qc.initialize(state3, [i for i in range(14,21)])
    qc.initialize(state4, [i for i in range(21,28)])


    qc.h(28)
    qc.barrier()
    for i in range(1,Nsite):
        for j in range(1,8):
            qc.cswap(28,j+(i-1)*7-1,j-1+21)
    qc.barrier()        
    qc.h(28)
    #qc.measure(28,0)
                 

    backend = Aer.get_backend('statevector_simulator')
    results = execute(qc,backend).result()
    s = results.get_statevector()
    a=s.probabilities([28])[0]
    b=s.probabilities([28])[1]
    return a-b
    


# In[4]:


def imcircuit():

    #from qiskit import QuantumCircuit, transpile, IBMQ
    #from qiskit import Aer, execute
    #from qiskit.visualization import plot_histogram
    #from qiskit.providers.aer import QasmSimulator
    #from qiskit.quantum_info import Statevector
    #import matplotlib as mpl

    Nsite = 4
    Nq = 29
    qc =  QuantumCircuit(Nq,1)

    state1 = np.zeros(2**7,dtype=complex)
    state2 = np.zeros(2**7,dtype=complex)
    state3 = np.zeros(2**7,dtype=complex)
    state4 = np.zeros(2**7,dtype=complex)



    state1[0:81] = matrix[0]
    state2[0:81] = matrix[1]
    state3[0:81] = matrix[2]
    state4[0:81] = matrix[3]


    qc.initialize(state1, [i for i in range(0,7)])
    qc.initialize(state2, [i for i in range(7,14)])
    qc.initialize(state3, [i for i in range(14,21)])
    qc.initialize(state4, [i for i in range(21,28)])


    qc.h(28)
    qc.barrier()
    for i in range(1,Nsite):
        for j in range(1,8):
            qc.cswap(28,j+(i-1)*7-1,j-1+21)
    qc.barrier()   
    qc.s(28)
    qc.h(28)
    #qc.measure(28,0)
                 

    backend = Aer.get_backend('statevector_simulator')
    results = execute(qc,backend).result()
    s = results.get_statevector()
    a=s.probabilities([28])[0]
    b=s.probabilities([28])[1]
    return b-a
    


    
    
    
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs


def eight_site_eig(sites,steps,start,stop):
    H=0
    Hsolo = 0
    Sx=(1/np.sqrt(2))*np.array([[0,1,0],[1,0,1],[0,1,0]])
    Sy=(1/np.sqrt(2))*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
    Sz=np.array([[1,0,0],[0,0,0],[0,0,-1]])
    I=np.array([[1,0,0],[0,1,0],[0,0,1]])
    Sup = np.sqrt(2)*np.array([[0,1,0],[0,0,1],[0,0,0]])
    Sdown = np.sqrt(2)*np.array([[0,0,0],[1,0,0],[0,1,0]])
    #nearest neighbor interaction XYZ
    output = np.array(np.zeros((steps,3**sites)),dtype=complex)

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
    berryphase = np.zeros(stop-start)
    index=0
    for j in range(start,stop):
        
        
        for i in range(steps):

            D = j/100
            phi = np.pi*i*2/(steps-1)
            #phi introduces the berry phase, for 4 sites, we have 0, 2pi/3, 4pi/3, 2pi
            H81 = 0.5*((np.exp(1j*phi)*minusterm) + (np.exp(-1j*phi)*plusterm))+zterm
            Ham=H+H81+D*Hsolo
            
            
            w,v = np.linalg.eig(Ham)
            #spa=csr_matrix(Ham)
            #w, v = eigs(spa,k=6,which='SR')
            
            wmin=np.min(w)
            ind=np.argmin(w)
            vmin=v[:,ind]
            output[i,:] = vmin
            
        result = 0
        
        for k in range(steps):
            if k != (steps-1):
                result = result+np.log(np.vdot(output[k],output[k+1]))
            elif k == (steps-1):
                result = result+np.log(np.vdot(output[k],output[0]))
    
        
        berryphase[index] = result.imag
        index=index+1
        
        
        
        
        
        
        
        
        
        
    return berryphase

    


datapoint = 80
x = np.linspace(0.4,1.2,datapoint)
berry = np.zeros(datapoint)
index = 0
for i in x:
    
    matrix = qberry(4,4,i)  #call the main function
    p = recircuit() + imcircuit()*1j
    berry[index] = cmath.phase(p)
    index=index+1
    


for i in range(len(berry)):
    if berry[i] < 0:
        berry[i] = -berry[i]

        
        

        

        
limit = np.linspace(0.4,1.2,80)



fig, ax1 = plt.subplots()


# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.55, 0.55, 0.3, 0.3]

ax2 = fig.add_axes([left, bottom, width, height])


a = eight_site_eig(4,4,40,120)
b = eight_site_eig(4,8,40,120)
c = eight_site_eig(4,12,40,120)
d = eight_site_eig(4,30,40,120)


for i in range(len(a)):
    if a[i]<0:
        a[i] = -a[i]
        
for i in range(len(b)):
    if b[i]<0:
        b[i] = -b[i]
        
for i in range(len(c)):
    if c[i]<0:
        c[i] = -c[i]

for i in range(len(d)):
    if d[i]<0:
        d[i] = -d[i]  
        

        
limit = np.linspace(0.4,1.2,80)

plt.rcParams["figure.figsize"] = (8,4)

ax1.plot(limit,berry,'o',label = 'k=4(quantum)')
ax2.plot(limit,a,color='red',label = 'k=4')
ax2.plot(limit,b,color='blue',label = 'k=8')
ax2.plot(limit,c,color='green',label = 'k=12')
ax2.plot(limit,d,color='black',label = 'k=30')



ax1.set_xlabel('D')
ax1.set_ylabel('Berry Phase')
ax2.set_xlabel('D')
ax2.set_ylabel('Berry Phase')


ax1.legend(loc='lower left')
ax2.legend(loc='best')        
        


plt.savefig('data')



# In[ ]:




