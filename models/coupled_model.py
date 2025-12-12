import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import  odeint, solve_ivp

from model_is import solve_imune, imune_parameters, imune_init_cond
from model_cv import solve_cv

# -----------------------------------------------------------------------------

def solve_modelo_acoplado_cv_im(params, cid=1, sid=0):
  """
  Args:
      params (_type_): _description_
      cid (int, optional): id do caso (sobrevivente ou nao). Defaults to 0.
      sid (int, optional): id da amostra. Defaults to 0.
  """
  
  # params['Cmax'] = valor
  
  Cmax = params['Cmax']

  # parametros e condicoes iniciais do imune
  pars = imune_parameters(caso=cid)
  y0_imune = imune_init_cond()

  # adiciona algumas condicoes iniciais nos parametros
  pars['Ap0']  = y0_imune[1] # Ap0
  pars['ThN0'] = y0_imune[4] # ThN0
  pars['TkN0'] = y0_imune[6] # TkN0
  pars['B0']   = y0_imune[8] # B0

  # parametros e condicoes iniciais do cardio
  y0cv = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

  # tempo total de simulacao do acoplamento
  tf = 120.0 # dias
  dt = 1.0   # dia 
  N = int(tf/dt)

  # janela de tempo para resolver o sistema cardiovascular
  dcv = 2.0 #dias
  tcv = np.arange(0.0, tf+dcv, dcv)
  sol = np.zeros(N)
  t = np.zeros(N)
  
  y1_anteriores = []
  yimune_anteriores = []

  tt = np.array(y0_imune)
  tt = tt.reshape(15,1)

  teste = np.hstack( ( tt, ) )

  # CV Janela 0
  pCovid = {}
  pCovid['krs'] = 1.0
  pCovid['krp'] = 1.0 
  pCovid['kt'] = 1.0
  pCovid['kea'] = 1.0
  
  print('Parametros covid:', end=' ')
  for key in pCovid:
    print(' %s : %.4f' % (key, pCovid[key]), end=',')
  print()

  print("Resolvendo o Modelo Cardiovascular passo: ", 0)
  y1cv = solve_cv(pCovid, 0.0, y0cv, step=0, sample_id=sid)
  y0cv = y1cv.copy()

  for i in range(len(tcv)-1):

    # -------------------------------------------------------------------------
    # Sistema Imune
    # -------------------------------------------------------------------------
    t_inicial = tcv[i] 
    t_final = tcv[i+1] 
    t_span = [t_inicial, t_final]

    #Janela de tempo para resolver o sistema imune
    tsi = np.arange(t_inicial, t_final + dt, dt) 
    nt = len(tsi) - 1
    
    print("Resolvendo o Modelo do Sistema Imune nos tempos: ", tsi)

    y1_imune = solve_imune(pars, t_span, y0_imune)    
    y0_imune = (y1_imune.y[:,-1]).copy()
    
    y1_anteriores.append(y1_imune.y[14,:])

    teste = np.hstack( (teste, y1_imune.y[:,:]) )

    citocina = y1_anteriores[-1]
    C = citocina[-1]
    print('  Citocina',C)
    
    # Calcula parametros do efeito do covid
    pCovid['krs'] = 1.0
    pCovid['krp'] = 1.0 + C/Cmax * 1.75 #resistencia pulmonar
    pCovid['kt'] = 1.0 - C/Cmax * 0.25 #frequencia cardiaca
    pCovid['kea'] = 1.0 - C/Cmax * 0.3 #elastancia
    print('  Parametros covid:', end=' ')
    for key in pCovid:
      print(' %s : %.4f' % (key, pCovid[key]), end=',')
    print()

    # -------------------------------------------------------------------------
    # Sistema CV
    # -------------------------------------------------------------------------
    print("Resolvendo o Modelo Cardiovascular passo: ", i+1)
    y1cv = solve_cv(pCovid, t_final, y0cv, step=i+1, sample_id=sid)
    y0cv = y1cv.copy()
   
  path = "outputs_%04d/" % sid
  if not os.path.exists(path):
      os.makedirs(path)

  nsize = np.shape(teste)[1]
  t = np.linspace(0.0,tf,nsize)
  dados_im = np.vstack( (t.reshape((1,nsize)), teste))
  
  # arq_out = "outputs_%04d/outputs_imune.txt" % sid  
  # np.savetxt(arq_out, dados_im)

  arq_out = "outputs_%04d/outputs_imune.npz" % sid
  np.savez(arq_out, data=dados_im)

  #return outputs

# -----------------------------------------------------------------------------

if __name__ == "__main__":

  Cmax = 12 #pg/mL
  p = {}
  p['Cmax'] = Cmax

  # resolve o modelo acoplado cardiovascular (cv) + imune (im)
  solve_modelo_acoplado_cv_im(p, cid=1, sid=0)
  

# end