import numpy as np
import math, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import pi, cos
from functions import *

# -----------------------------------------------------------------------------

def ElastanceAtrium(t,EaM,Eam,Tar,tac,Tac,T):
    """
    Input the following variables:
    t = time from 0 to T
    T = length of cardiac cycle
    """
    if t<=Tar:
        Ea = (EaM-Eam)*(1-cos(pi*(t-Tar)/(T-Tac+Tar)))/2+Eam
    elif t <= tac:
        Ea = Eam
    elif t <= Tac:
        Ea = (EaM-Eam)*(1-cos(pi*(t-tac)/(Tac-tac)))/2+Eam
    else:
        Ea = (EaM-Eam)*(1+cos(pi*(t-Tac)/(T-Tac+Tar)))/2+Eam
    return Ea

# -----------------------------------------------------------------------------    

def ElastanceVentricle(t,EvM,Evm,Tvc,Tvr):
    """
    Input the following variables:
    t = time from 0 to T
    T = heart rate
    """
    if t<=Tvc:
        Elv = (EvM-Evm)*(1-cos(pi*t/Tvc))/2 + Evm
    elif t <= Tvr:
        Elv = (EvM-Evm)*(1+cos(pi*(t-Tvc)/(Tvr-Tvc)))/2 + Evm
    else:
        Elv = Evm  
    return Elv

# -----------------------------------------------------------------------------

def Modelo(y,t,pars,KEa,KRs,KRp,KT):

    Vsa = y[0]  # systemic artery volume
    Vsv = y[1]  # systemic venous volume
    Vpa = y[2]  # pulmonary artery volume
    Vpv = y[3]  # pulmonary venous volume
    Vra = y[4]  # right atria volume
    Vla = y[5]  # left atria volume
    Vrv = y[6]  # right ventricular volume
    Vlv = y[7]  # left ventricular volume

    # Parameters --------------------------------

    # Resistances
    Rs   = pars['Rs']   # systemic arteries
    Rp   = pars['Rp']   # pulmonary arteries
    Rava = pars['Rava'] # aortic valve
    Rmva = pars['Rmva'] # mitral valve
    Rpva = pars['Rpva'] # pulmonary valve
    Rtva = pars['Rtva'] # tricuspid valve
    Rpv  = pars['Rpv']  # pulmonary veins
    Rsv  = pars['Rsv']  # systemic veins

    # Compliances
    Csa = pars['Csa']   # systemic arteries
    Csv = pars['Csv']   # systemic veins
    Cpa = pars['Cpa']   # pulmonary arteries
    Cpv = pars['Cpv']   # pulmonary veins
    
    # Heart Elastances
    EMra = pars['EMra'] # right atrium max 
    Emra = pars['Emra'] # right atrium min 
    EMla = pars['EMla'] # left atrium max 
    Emla = pars['Emla'] # left atrium min 
    EMrv = pars['EMrv'] # right ventricle max 
    Emrv = pars['Emrv'] # right ventricle min 
    EMlv = pars['EMlv'] # left ventricle max
    Emlv = pars['Emlv'] # left ventricle min 

    # Timings
    Trra = pars['Trra'] # end right atrium relaxation
    tcra = pars['tcra'] # right atrium begin contraction
    Tcra = pars['Tcra'] # right atrium end contraction
    Trla = pars['Trla'] # end left atrium relaxation
    tcla = pars['tcla'] # left atrium begin contraction
    Tcla = pars['Tcla'] # left atrium end contraction
    Tcrv = pars['Tcrv'] # right ventricle contraction
    Trrv = pars['Trrv'] # right ventricle relaxation
    Tclv = pars['Tclv'] # left ventricle contraction
    Trlv = pars['Trlv'] # left ventricle relaxation

    T = pars['T']
    tst = t - (t % T)

    # Timevarying elastance
    Era = ElastanceAtrium(t-tst,EMra,Emra,Trra,tcra,Tcra,T) # Elastance right atrium
    Ela = ElastanceAtrium(t-tst,EMla,Emla,Trla,tcla,Tcla,T) # Elastance left atrium
    Erv = ElastanceVentricle(t-tst,EMrv,Emrv,Tcrv,Trrv)     # Elastance right ventricle
    Elv = ElastanceVentricle(t-tst,EMlv,Emlv,Tclv,Trlv)     # Elastance left ventricle

    psa = Vsa/Csa # systemic arteries 
    psv = Vsv/Csv # systemic veins 
    ppa = Vpa/Cpa # pressure pulmonary artery
    ppv = Vpv/Cpv # pressure pulmonary vein

    # Pressure 
    pra  = Era*Vra # pressure right atria
    prv  = Erv*Vrv # pressure right ventricle
    pla  = Ela*Vla # pressure left atria
    plv  = Elv*Vlv # pressure left ventricle

    if psv > pra:
        qsv = (psv-pra)/Rsv
    else:
        qsv = 0

    # Flow through valves
    
    # flow through tricuspid valve
    if pra > prv:  
        qtva = (pra-prv)/Rtva  # valve open
    else:
        qtva = 0              # valve closed

    # flow through mitral valve
    if pla > plv: 
        qmva = (pla-plv)/Rmva # valve open
    else:
        qmva = 0              # valve closed

    # flow through aortic valve 
    if plv > psa: 
        #dqava = (plv-psa - qava*Rava)/Lava; # valve open
        qava = (plv-psa)/Rava
    else:
    #   dqava = 0;                           # valve closed
        qava  = 0

    # flow through pulmonary valve
    if prv > ppa: 
    #%   dqpva = (prv-ppa-qpva*Rpva)/Lpva; % valve open
        qpva = (prv-ppa)/Rpva
    else:
    #%   dqpva = 0;                     % valve closed
        qpva  = 0

    # Peripheral flows
    qs = (psa - psv) / Rs # Systemic periphery
    qp = (ppa - ppv) / Rp # Pulmonary periphery

    # Venous flows (with venous valves - flow cannot go backwards into the veins)
    if ppv > pla:
        qpv = (ppv-pla)/Rpv
    else:
        qpv = 0

    if psv > pra:
        qsv = (psv-pra)/Rsv
    else:
        qsv = 0

    # Differential Equations
    dVsa = qava - qs
    dVsv = qs   - qsv
    dVpa = qpva - qp
    dVpv = qp   - qpv
    dVra = qsv  - qtva 
    dVla = qpv  - qmva
    dVrv = qtva - qpva
    dVlv = qmva - qava

    dydt = [dVsa, dVsv, dVpa, dVpv, dVra, dVla, dVrv, dVlv]

    return dydt

# -----------------------------------------------------------------------------

def ReadParsInit(filename,kRs,kRp,kEa,kT):

    # Parametros
    # [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
    #  Csa Csv Cpa Cpv ...                    % 9-12
    #  EMra Emra EMla Emla ...                % 13-16
    #  EMrv Emrv EMlv Emlv ...                % 17-20
    #  Trra tcra Tcra ...                     % 21-23
    #  Tcrv Trrv]';                           % 24-25

    #filename = 'paciente1.xlsx'  
    #filename = 'paciente-normotensive.xlsx'

    dados = pd.read_excel(filename)  
    
    pars = dados['pars'] # parameteres
    init = dados['init'] # initial conditions

    # dictionary containing the parameters
    p = {} 

    p['Rs'] = kRs * pars[0]
    p['Rp'] = kRp * pars[1]
    p['Rava'] = pars[2]
    p['Rmva'] = pars[3]
    p['Rpva'] = pars[4]
    p['Rtva']  = pars[5]
    p['Rpv'] = kRp * pars[6]
    p['Rsv']  = kRs * pars[7]

    p['Csa'] = pars[8]
    p['Csv'] = pars[9]
    p['Cpa'] = pars[10]
    p['Cpv'] = pars[11]

    p['EMra'] = kEa * pars[12]
    p['Emra'] = kEa * pars[13]
    p['EMla'] = kEa * pars[14]
    p['Emla'] = kEa * pars[15]

    p['EMrv'] = kEa * pars[16]
    p['Emrv'] = kEa * pars[17]
    p['EMlv'] = kEa * pars[18]
    p['Emlv'] = kEa * pars[19]

    Trra = kT * pars[20]
    tcra = kT * pars[21]
    Tcra = kT * pars[22]
    Tcrv = kT * pars[23]
    Trrv = kT * pars[24]

    p['T'] = 1.0*kT

    # Timing parameters
    p['Trra'] = Trra                   # end right atrium relaxation
    p['tcra'] = p['Trra'] + tcra       # right atrium begin contraction
    p['Tcra'] = p['tcra'] + Tcra       # right atrium end contraction

    p['Trla'] = p['Trra']*(1.01)       # end left atrium relaxation
    p['tcla'] = p['tcra']*(1.05)       # left atrium begin contraction
    p['Tcla'] = p['Tcra']              # left atrium end contraction

    p['Tcrv'] = Tcrv                   # right ventricle contraction
    p['Trrv'] = p['Tcrv'] + Trrv       # right ventricle relaxation
    
    p['Tclv'] = p['Tcrv'] *(0.95)      # left ventricle contraction
    p['Trlv'] = p['Trrv']              # left ventricle relaxation
 
    # initial conditions
    Vsa0 = init[0]  # systemic artery volume
    Vsv0 = init[1]  # systemic venous volume
    Vpa0 = init[2]  # pulmonary artery volume
    Vpv0 = init[3]  # pulmonary venous volume
    Vra0 = init[4]  # right atria volume
    Vla0 = init[5]  # left atria volume
    Vrv0 = init[6]  # right ventricular volume
    Vlv0 = init[7]  # left ventricular volume
    y0 = [Vsa0, Vsv0, Vpa0, Vpv0, Vra0, Vla0, Vrv0, Vlv0]

    return p, y0

#-------------------------------------------------------------------------------

def solve_cv(p, t, y0, step = 0, sample_id = 0):

    # normal: kT = kEa = kRs = kRp = 1.0
    # ENMC: kT = 0.75, kEa = 0.70, kRs = 1.0, kRp = 2.75

    filename = 'paciente-normotensive.xlsx'
    #path = 'outputs/'

    kT  = p['kt']   # heart rate changes
    kEa = p['kea']  # elastance changes
    kRs = p['krs']  # systemic resistance changes
    kRp = p['krp']  # pulmonary resistance changes

    pars, y0 = ReadParsInit(filename,kRs,kRp,kEa,kT)
    # print(kRs,kRp,kEa,kT)    
    # print(pars)

    # y0 -> sempre pega o mesmo (dos dados do arquivo)

    # tempo
    nc = 30
    ns = 1000
    T = pars['T']
    N = nc*ns
    tf = nc*T
    t = np.linspace(0, tf, nc*ns)
    k = (nc-1)*ns # indice do ultimo ciclo

    # print(y0)

    # solve the model
    solution = odeint(Modelo, y0, t, args=(pars,kEa,kRs,kRp,kT))

    # get solution
    Vsa = solution[:,0]
    Vsv = solution[:,1]
    Vpa = solution[:,2]
    Vpv = solution[:,3]
    Vra = solution[:,4]
    Vla = solution[:,5]
    Vrv = solution[:,6]
    Vlv = solution[:,7]

    # -------------------------------------------------------------------------

    # post-processing
    Era = np.zeros(N)
    Ela = np.zeros(N)
    Erv = np.zeros(N)
    Elv = np.zeros(N)

    psa = Vsa/pars['Csa'] # systemic arteries
    psv = Vsv/pars['Csv'] # systemic veins
    ppa = Vpa/pars['Cpa'] # pulmonary arteries
    ppv = Vpv/pars['Cpv'] # pulmonary veins
    
    # Heart elastance
    for i in range(N):
        tst = t[i] - (t[i] % T)
        Era[i] = ElastanceAtrium(t[i]-tst,pars['EMra'],pars['Emra'],pars['Trra'],pars['tcra'],pars['Tcra'],pars['T']) # Right atrium elastance
        Ela[i] = ElastanceAtrium(t[i]-tst,pars['EMla'],pars['Emla'],pars['Trla'],pars['tcla'],pars['Tcla'],pars['T']) # Left atrium elastance        
        Erv[i] = ElastanceVentricle(t[i]-tst,pars['EMrv'],pars['Emrv'],pars['Tcrv'],pars['Trrv'])     # Right ventricle elastance
        Elv[i] = ElastanceVentricle(t[i]-tst,pars['EMlv'],pars['Emlv'],pars['Tclv'],pars['Trlv'])     # Right ventricle elastance
    
    pra = Era * Vra # Right atrium pressure
    pla = Ela * Vla # Left atrium pressure
    prv = Erv * Vrv # Right ventricle pressure
    plv = Elv * Vlv # Left ventricle pressure

    #Salvamento de dados
    path = "outputs_%04d/" % sample_id
    if not os.path.exists(path):
        os.makedirs(path)

    #np.savetxt(path + '/outputs.txt', 
        #np.transpose([t, Vsa, Vsv, Vpa, Vpv,Vra,Vla,Vrv,Vlv,
                         #psa,psv,ppa,ppv,pra,pla,prv,plv, 
                         #Era,Ela,Erv,Elv]))

    dados_saida = np.transpose([t[k:], Vsa[k:], Vsv[k:], Vpa[k:], Vpv[k:],Vra[k:],Vla[k:],Vrv[k:],Vlv[k:],
                             psa[k:],psv[k:],ppa[k:],ppv[k:],pra[k:],pla[k:],prv[k:],plv[k:], 
                             Era[k:],Ela[k:],Erv[k:],Elv[k:]])

    # arq_out = 'outputs_%04d/outputs_cv_%02d.txt' % (sample_id, step)   
    # np.savetxt(arq_out, dados_saida)

    arq_out = 'outputs_%04d/outputs_cv_%02d.npz' % (sample_id, step)
    np.savez_compressed(arq_out, data=dados_saida)
    
    CalcOutputs(k,t,T,psa,psv,ppa,ppv,Vlv,plv,Vrv,prv,step,kT)
    ImprimeCalcOutputs(step)
    #PlotAll(t[k:],plv[k:],prv[k:],pla[k:],pra[k:],psa[k:],psv[k:],ppa[k:],ppv[k:],path)
    #PlotElastances(t[k:],Elv[k:],Ela[k:],Erv[k:],Era[k:],path)
    #PlotPVLoops(Vlv[k:],plv[k:],Vrv[k:],prv[k:],step)  
    #PlotVolumes(k,t,solution,path)
    #PlotPressures(t[k:],plv[k:],psa[k:],psv[k:],ppa[k:],ppv[k:],path)

    return np.array([Vsa[-1], Vsv[-1], Vpa[-1], Vpv[-1], Vra[-1], Vla[-1], Vrv[-1], Vlv[-1]])

#-------------------------------------------------------------------------------

'''if __name__ == "__main__":

    # parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='file', type=str, help='input filename')
    parser.add_argument('-o', dest='outdir', type=str, help='output directory')
    parser.add_argument('-kt', dest='kt',  type=float, default=1.0, help='factor to change T')
    parser.add_argument('-kea', dest='kea', type=float, default=1.0, help='factor to change Ea')
    parser.add_argument('-krs', dest='krs', type=float, default=1.0, help='factor to change Rs')
    parser.add_argument('-krp', dest='krp', type=float, default=1.0, help='factor to change Rp')
    args = parser.parse_args()

    # normal: kT = kEa = kRs = kRp = 1.0
    # ENMC: kT = 0.75, kEa = 0.70, kRs = 1.0, kRp = 2.75

    kT  = args.kt   # heart rate changes
    kEa = args.kea  # elastance changes
    kRs = args.krs  # systemic resistance changes
    kRp = args.krp  # pulmonary resistance changes
    
    filename = 'paciente-normotensive.xlsx'
    if(args.file is not None):
        filename = args.file

    #path = 'outputs/'
    #if(args.outdir is not None):
        #path = args.outdir

    # read the parameters and initial conditions from XLS file
    pars, y0 = ReadParsInit(filename,kRs,kRp,kEa,kT)
        
    resolve_cv(pars, 0, y0)'''



# end