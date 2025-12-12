import os, sys
sys.path.append('models/')

import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import  odeint, solve_ivp

from coupled_model import * 
from model_is import solve_imune
from model_cv import solve_cv

if __name__ == "__main__":

    nsa  = 1000 #1000 # numero de amostras
    ninp = 1   # numero de parametros
    nout = 8   # numero de saidas (qois)

    caso = 0 # survivors
    # caso = 1 # non-survivors

    # Caso 1
    #Cmax = cp.Normal(12,1.2) # cov = std/med = 10% = 0.10

    # Caso 2
    # ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8531346/pdf/41598_2021_Article_190.pdf
    Cmax = cp.Uniform(8.07, 26.98) 

    #distribution = cp.J(d_kEa,d_kRp,d_kT,Cmax)
    distribution = cp.J(Cmax)
    
    samples = distribution.sample(nsa)
    evals = np.empty((nsa,nout))
    y0cv = np.empty((nsa,nout))
    
    for i in range(nsa):

        print('\nCaso %d\n' % i)
        r = {}
        r['Cmax'] = samples[i]
        # r['Cmax'] = 12.0
        #evals[i,:] = solve_modelo_acoplado_cv_im(r, sid=0)
        evals[i,:] = solve_modelo_acoplado_cv_im(r, cid=caso, sid=i)