import numpy as np
import math, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
#plt.style.use('science')
from scipy.integrate import odeint
from math import pi, cos

from model_cv import *

# -----------------------------------------------------------------------------
def CalcOutputs(k,t,T,psa,psv,ppa,ppv,Vlv,plv,Vrv,prv,id,kt):

    # Systemic Circulation
    
    #Systemic Artery
    Psa_min = psa[k:].min()
    Psa_max = psa[k:].max()
    Psa_mean = (1.0/T)*np.trapz(psa[k:],t[k:])

    #Systemic Vein
    Psv_min = psv[k:].min()
    Psv_max = psv[k:].max()
    Psv_mean = (1.0/T)*np.trapz(psv[k:],t[k:])   

    #Pulmonary Circulation

    #Pulmonary Artery
    Ppa_min = ppa[k:].min()
    Ppa_max = ppa[k:].max()
    Ppa_mean = (1.0/T)*np.trapz(ppa[k:],t[k:])

    #Pulmonary Vein
    Ppv_min = ppv[k:].min()
    Ppv_max = ppv[k:].max()
    Ppv_mean = (1.0/T)*np.trapz(ppv[k:],t[k:])

    Plv_min = plv[k:].min()
    Plv_max = plv[k:].max()
    Plv_mean = (1.0/T)*np.trapz(plv[k:],t[k:])

    #Left Ventricle
    LV_EDV = Vlv[k:].max()
    LV_ESV = Vlv[k:].min()
    LV_SV = LV_EDV-LV_ESV
    LV_EF = (LV_SV/LV_EDV)*100
    LV_Pmax = plv[k:].max()
    LV_CO = (LV_SV*(60.0/T*kt))

    # Right Ventricule
    RV_EDV = Vrv[k:].max()
    RV_ESV = Vrv[k:].min()
    RV_SV = RV_EDV - RV_ESV
    RV_EF = (RV_SV/RV_EDV)*100
    RV_Pmax = prv[k:].max()
    RV_CO = (RV_SV*(60.0/T*kt))

    #Salvamento de dados
    path = 'CalcOutputs/'

    if not os.path.exists(path):
        os.makedirs(path)
   
    arq =  path + 'outputs' + str(id) + '.txt'
    
    hdr = "Psa_min, Psa_max, Psa_mean, Psv_min, Psv_max, Psv_mean, Ppa_min, Ppa_max, Ppa_mean, Ppv_min, Ppv_max, Ppv_mean, Plv_min, Plv_max, Plv_mean, LV_EDV, LV_ESV, LV_SV, LV_EF, LV_Pmax, LV_CO, RV_EDV, RV_ESV, RV_SV, RV_EF, RV_Pmax, RV_CO"
    out = np.array([Psa_min, Psa_max, Psa_mean, Psv_min, Psv_max, Psv_mean, Ppa_min, Ppa_max, Ppa_mean, Ppv_min, Ppv_max, Ppv_mean, Plv_min, Plv_max, Plv_mean, LV_EDV, LV_ESV, LV_SV, LV_EF, LV_Pmax, LV_CO, RV_EDV, RV_ESV, RV_SV, RV_EF, RV_Pmax, RV_CO])
    np.savetxt(arq, out, header=hdr)

# -----------------------------------------------------------------------------

def ImprimeCalcOutputs(id):
    arq = 'CalcOutputs/outputs' + str(id) + '.txt'
    dados = np.loadtxt(arq)
    
    Q = {}
    #Systemic Artery
    Q['Psa_min'] = dados[0]
    Q['Psa_max'] = dados[1]
    Q['Psa_mean'] = dados[2]

    #Systemic Vein
    Q['Psv_min'] = dados[3]
    Q['Psv_max'] = dados[4]
    Q['Psv_mean'] = dados[5] 

    #Pulmonary Circulation

    #Pulmonary Artery
    Q['Ppa_min'] = dados[6]
    Q['Ppa_max'] = dados[7]
    Q['Ppa_mean'] = dados[8]

    #Pulmonary Vein
    Q['Ppv_min'] = dados[9]
    Q['Ppv_max'] = dados[10]
    Q['Ppv_mean'] = dados[11]

    Q['Plv_min'] = dados[12]
    Q['Plv_max'] = dados[13]
    Q['Plv_mean'] = dados[14]

    #Left Ventricle
    Q['LV_EDV'] = dados[15]
    Q['LV_ESV'] = dados[16]
    Q['LV_SV'] = dados[17]
    Q['LV_EF'] = dados[18]
    Q['LV_Pmax'] = dados[19]
    Q['LV_CO'] = dados[20]

    # Right Ventricule
    Q['RV_EDV'] = dados[21]
    Q['RV_ESV'] = dados[22]
    Q['RV_SV'] = dados[23]
    Q['RV_EF'] = dados[24]
    Q['RV_Pmax'] = dados[25]
    Q['RV_CO'] = dados[26]

    #print(' ', Q)

# -----------------------------------------------------------------------------

def PlotAll(t,plv,prv,pla,pra,psa,psv,ppa,ppv,path):
    
    fig, axs = plt.subplots(3, 3)

    # linha 1
    axs[0, 0].plot(t, pla)
    axs[0, 0].set(xlabel='time [s]', ylabel='pla')

    axs[0, 1].plot(t, plv, 'tab:orange')
    axs[0, 1].set(xlabel='time [s]', ylabel='plv')

    axs[0, 2].plot(t, psa, 'tab:orange')
    axs[0, 2].set(xlabel='time [s]', ylabel='psa')
    
    # linha 2
    axs[1, 0].plot(t, psv, 'tab:green')
    axs[1, 0].set(xlabel='time [s]', ylabel='psv')
    
    axs[1, 1].plot(t, pra, 'tab:red')
    axs[1, 1].set(xlabel='time [s]', ylabel='pra')

    axs[1, 2].plot(t, prv, 'tab:red')
    axs[1, 2].set(xlabel='time [s]', ylabel='prv')

    # linha 3
    axs[2, 0].plot(t, ppa, 'tab:green')
    axs[2, 0].set(xlabel='time [s]', ylabel='ppa')
    
    axs[2, 1].plot(t, ppv, 'tab:red')
    axs[2, 1].set(xlabel='time [s]', ylabel='ppv')

    #axs[2, 1].plot(t, psa, 'tab:red')
    #axs[2, 2].set_title('Axis [2, 2]')

    plt.tight_layout()
    plt.savefig(path+'fig_pressures.pdf', format='pdf', dpi=300)
    plt.savefig(path+'fig_pressures.png', format='png', dpi=300)
    plt.show()

# -----------------------------------------------------------------------------

def PlotElastances(t,Elv,Ela,Erv,Era,path):
    plt.plot(t, Elv, label ='Elv')
    plt.plot(t, Ela, label ='Ela')
    plt.plot(t, Erv, label ='Erv')
    plt.plot(t, Era, label ='Era')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('elastance')
    plt.tight_layout()
    plt.savefig(path+'fig_elastance.png', format='png', dpi=300)
    plt.show()

# -----------------------------------------------------------------------------

def PlotPressures(t,plv,psa,psv,ppa,ppv,path):
    
    plt.plot(t,plv, label ='PLV')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('PLV [mmHg]')
    plt.tight_layout()
    plt.savefig(path+'fig_plv.png', format='png', dpi=300)
    plt.show()

    plt.plot(t, psa, label='Psa') 
    #plt.legend(loc='best')
    plt.xlabel('tempo [s]')
    plt.ylabel('pressão [mmHg]')
    plt.title(r'$P_{sa}$')
    plt.tight_layout()
    plt.savefig(path+'fig_psa.pdf', format='pdf', dpi=300)
    plt.show()

    plt.plot(t, psv, label='Psv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('pressure [mmHg]')
    plt.savefig(path+'fig_psv.png', format='png', dpi=300)
    plt.show()

    plt.plot(t, ppa, label='Ppa')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('pressure [mmHg]')
    plt.savefig(path+'fig_ppa.png', format='png', dpi=300)
    plt.show()

    plt.plot(t, ppv, label='Ppv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('pressure [mmHg]')
    plt.savefig(path+'fig_ppv.png', format='png', dpi=300)
    plt.show()

# -----------------------------------------------------------------------------

def PlotVolumes(k,t,solution,path):
        
    plt.plot(t[k:], solution[k:,0], label='Vsa')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vsa.png', format='png', dpi=300)

    plt.plot(t[k:], solution[k:,1], label='Vsv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vsv.png', format='png', dpi=300)

    plt.plot(t[k:], solution[k:,2], label='Vpa')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vpa.png', format='png', dpi=300)

    plt.plot(t[k:], solution[k:,3], label='Vpv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vpv.png', format='png', dpi=300)
    plt.show()

    plt.plot(t[k:], solution[k:,4], label='Vra')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vra.png', format='png', dpi=300)
    plt.show()

    plt.plot(t[k:], solution[k:,5], label='Vla')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vla.png', format='png', dpi=300)
    plt.show()

    plt.plot(t[k:], solution[k:,6], label='Vrv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vrv.png', format='png', dpi=300)
    plt.show()

    plt.plot(t[k:], solution[k:,7], label='Vlv')
    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('volume [mL]')
    plt.savefig(path+'fig_vlv.png', format='png', dpi=300)
    plt.show()

# -----------------------------------------------------------------------------

def PlotPVLoops(Vlv,plv,Vrv,prv,id):

    #plt.plot(Vrv, prv)
    #plt.xlabel('Vrv [mL]')
    #plt.ylabel('Prv [mmHg]')
    #plt.tight_layout()
    #plt.savefig(path+'fig_pvloop_rv.png', format='png', dpi=300)
    #plt.show()

    #plt.plot(Vlv, plv)
    #plt.xlabel('Vlv [mL]')
    #plt.ylabel('Plv[mmHg]')
    #plt.tight_layout()
    #plt.savefig(path+'fig_pvloop_lv.png', format='png', dpi=300)
    #plt.show()

    #Salvamento de dados
    path = 'PlotPVLoops/'

    if not os.path.exists(path):
        os.makedirs(path)
   
    arq =  path + 'outputs' + str(id) + '.txt'
    
    hdr = "Vlv,plv,Vrv,prv"
    out = np.array([Vlv,plv,Vrv,prv])
    np.savetxt(arq, out, header=hdr)