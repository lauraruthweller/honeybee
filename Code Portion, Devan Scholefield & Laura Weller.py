#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Submitted by Devan Scholefield and Laura Weller


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

L=2000    
W=np.array([4000,10000,27000])
N=np.arange(0,40000,100)

for w_value in W:
    for N_value in N:
        Y=(L*N_value)/(N_value+w_value)
        plt.plot(N_value,Y,'.', color='black')


# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def J_fun(m):
    L = 2000
    α = 0.25
    σ = 0.75
    w = 27000
    J = 0.5*(((α/m)-(σ/m)-1)+(np.sqrt(((((α/m)-(σ/m)-1))**2)+4*(α/m))))
    return J


def Fss(m,J):
    L = 2000
    w = 27000
    F0 = ((L/m)-(w*(J/(J+1))))
    return F0


def Hss(J,F0):
    H0 = ((1/J)*(F0))
    return H0

#Solving for m so that F and H are 0
m=0.355196688
print(J_fun(m))
print(Fss(m,J_fun(m)))
print(Hss(J_fun(m),Fss(m,J_fun(m))))



# In[4]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def ode(state,t):
    H = state[0]
    F = state[1]
    L = 2000
    w = 27000
    E = L*((H + F)/(w + H + F))
    α = 0.25
    σ = 0.75
    R = α-σ*(F/(H+F))
    
    
    dH = E-(H*(R))
    dF = (H*R)-(m*F)
    
    rstate=[dH,dF]
    
    return rstate


# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

tmax=100
dt=1
t=np.arange(0,tmax,dt)
m=0.24

x=Hss(J_fun(m),Fss(m,J_fun(m)))
y=Fss(m,J_fun(m))
plt.scatter(x,y,c='red')

xvalue=np.arange(start=0,stop=18000,step=3000)
yvalue=np.arange(start=0,stop=7200,step=1200)

for xvar in xvalue:
    for yvar in yvalue:
        state0=[xvar, yvar]
        state=odeint(ode,state0,t)
        x=state[:,0]
        y=state[:,1]
        plt.plot(x,y,'-',c='darkgrey')
plt.xlabel("H")
plt.ylabel("F")


# In[6]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

tmax=100
dt=1 
t=np.arange(0,tmax,dt)
m=0.4


plt.scatter(0,0,c='red')


xvalue=np.arange(start=0,stop=18000,step=3000)
yvalue=np.arange(start=0,stop=7200,step=1200)

for xvar in xvalue:
    for yvar in yvalue:
        state0=[xvar, yvar]
        state=odeint(ode,state0,t)
        x=state[:,0]
        y=state[:,1]
        plt.plot(x,y,'-',c='darkgrey')
plt.xlabel("H")
plt.ylabel("F")


# In[7]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def ode1(state,t):
    H = state[0]
    F = state[1]
    L = 2000
    m = 0.4
    w = 27000
    E = L*((H + F)/(w + H + F))
    α = 0.25
    σ = 0.75
    R = α-σ*(F/(H+F))
    
    dH = E-(H*(R))
    dF = (H*R)-(m*F)
    
    rstate=[dH,dF]
    
    return rstate



def ode2(state,t):
    H = state[0]
    F = state[1]
    λ = 0.059
    m1 = 0.6
    L = 2000
    w = 27000
    m = 0.6
    E = L*((H + F)/(w + H + F))
    α = 0.25
    σ = 0.75
    R = α-σ*(F/(H+F))
    dH = E-(H*(R))
    if R < 0:
        M = (m1*(R**2))/(λ+(R**2))
        dF = (H*R)-(M*F)
    else:
        dF = (H*R)-(m*F)
    
    rstate=[dH,dF]
    
    return rstate

tmax=200
dt=1 # time step
t=np.arange(0,tmax,dt)


state0=[14000,7000] 
state=odeint(ode1,state0,t)
state1=odeint(ode2,state0,t)
x=state[:,0]
y=state[:,1]
x1=state1[:,0]
y1=state1[:,1]
plt.plot(t,(x+y),':',c='blue')
plt.plot(t,(x1+y1),'-',c='red')
plt.xlabel('time (days)')
plt.ylabel('Number of bees in colony (N)')


# In[ ]:




