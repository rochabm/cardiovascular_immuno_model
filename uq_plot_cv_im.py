import os, sys
sys.path.append('models/')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from functions import *
import scienceplots
plt.style.use('nature')

# -----------------------------------------------------------------------------

def CalcQoIs(vlv_mat,plv_mat,vrv_mat,prv_mat,psa_mat,ppa_mat,pra_mat,sid,nc):

    print('qois from sample %d' % sid)
    qoi_data = {}
    qoi_data['LV_EDV']  = np.zeros(nc)
    qoi_data['LV_ESV']  = np.zeros(nc)
    qoi_data['LV_SV']   = np.zeros(nc)
    qoi_data['LV_EF']   = np.zeros(nc)
    qoi_data['LV_Pmax'] = np.zeros(nc)
    qoi_data['LV_Pmin'] = np.zeros(nc)

    qoi_data['RV_EDV']  = np.zeros(nc)
    qoi_data['RV_ESV']  = np.zeros(nc)
    qoi_data['RV_SV']   = np.zeros(nc)
    qoi_data['RV_EF']   = np.zeros(nc)
    qoi_data['RV_Pmax'] = np.zeros(nc)
    qoi_data['RV_Pmin'] = np.zeros(nc)

    # artéria sistêmica
    qoi_data['Psa_min'] = np.zeros(nc)
    qoi_data['Psa_max'] = np.zeros(nc)

    # artéria pulmonar
    qoi_data['Pmin'] = np.zeros(nc)
    qoi_data['Pmax'] = np.zeros(nc)
    
    # atrio direito
    qoi_data['Pra_min'] = np.zeros(nc)
    qoi_data['Pra_max'] = np.zeros(nc)

    for j in range(nc):

        # criar vetores para armazenar as quantidades fisiológicas
        print(f"Quantidades fisiológicas do passo {j}")

        # left ventricle
        qoi_data['LV_EDV'][j]  = vlv_mat[sid, j, :].max()
        qoi_data['LV_ESV'][j]  = vlv_mat[sid, j, :].min()
        qoi_data['LV_SV'][j]   = qoi_data['LV_EDV'][j] - qoi_data['LV_ESV'][j]
        qoi_data['LV_EF'][j]   = (qoi_data['LV_SV'][j]/qoi_data['LV_EDV'][j]) * 100
        qoi_data['LV_Pmax'][j] = plv_mat[sid, j, :].max()
        qoi_data['LV_Pmin'][j] = plv_mat[sid, j, :].min()

        print( qoi_data['LV_EF'][j]  )

        # right ventricule
        qoi_data['RV_EDV'][j]  = vrv_mat[sid, j, :].max()
        qoi_data['RV_ESV'][j]  = vrv_mat[sid, j, :].min()
        qoi_data['RV_SV'][j]   = qoi_data['RV_EDV'][j] - qoi_data['RV_ESV'][j]
        qoi_data['RV_EF'][j]   = (qoi_data['RV_SV'][j]/qoi_data['RV_EDV'][j]) * 100
        qoi_data['RV_Pmax'][j] = prv_mat[sid, j, :].max()
        qoi_data['RV_Pmin'][j] = prv_mat[sid, j, :].min()

        print(' LV_EDV   = %.2f' % qoi_data['LV_EDV'][j])
        print(' LV_ESV   = %.2f' % qoi_data['LV_ESV'][j])
        print(' LV_SV   = %.2f'  % qoi_data['LV_SV'][j])
        print(' LV_EF    = %.2f' % qoi_data['LV_EF'][j])
        print(' LV_Pmax  = %.2f' % qoi_data['LV_Pmax'][j])
        print(' LV_Pmin  = %.2f' % qoi_data['LV_Pmin'][j])
        print(' RV_EDV   = %.2f' % qoi_data['RV_EDV'][j])
        print(' RV_ESV   = %.2f' % qoi_data['RV_ESV'][j])
        print(' RV_SV   = %.2f'  % qoi_data['RV_SV'][j])
        print(' RV_EF    = %.2f' % qoi_data['RV_EF'][j])
        print(' RV_Pmax  = %.2f' % qoi_data['RV_Pmax'][j])
        print(' RV_Pmin  = %.2f' % qoi_data['RV_Pmin'][j])

        # systemic circulation
        qoi_data['Psa_min'][j] = psa_mat[sid, j, :].min()
        qoi_data['Psa_max'][j] = psa_mat[sid, j, :].max()
        print(' Psa_min  = %.2f' % qoi_data['Psa_min'][j])
        print(' Psa_max  = %.2f' % qoi_data['Psa_max'][j])

        # pulmonary circulation
        qoi_data['Pmin'][j] = ppa_mat[sid, j, :].min()
        qoi_data['Pmax'][j] = ppa_mat[sid, j, :].max()
        print(' PAP_min = %.2f' % qoi_data['Pmin'][j])
        print(' PAP_max = %.2f' % qoi_data['Pmax'][j])

        # right atrium
        qoi_data['Pra_min'][j] = pra_mat[sid, j, :].min()
        qoi_data['Pra_max'][j] = pra_mat[sid, j, :].max()
        print(' Pra_min  = %.2f' % qoi_data['Pra_min'][j])
        print(' Pra_max  = %.2f' % qoi_data['Pra_max'][j])

    #return LV_EDV, LV_ESV, LV_SV, LV_EF, LV_Pmax, RV_EDV, RV_ESV, RV_SV, RV_EF, RV_Pmax, RV_Pmin, Psa_min, Psa_max, Pmin, Pmax, Pra_min, Pra_max
    return qoi_data

# -----------------------------------------------------------------------------

def plot_pvloop(v_mat, p_mat, label, ciclo, cor='k'):
    pl = 5
    pm = 50
    pu = 95
    vrv_p5  = np.percentile(v_mat[:, ciclo, :], pl, axis=0)
    prv_p5  = np.percentile(p_mat[:, ciclo, :], pl, axis=0)    
    vrv_p50 = np.percentile(v_mat[:, ciclo, :], pm, axis=0)
    prv_p50 = np.percentile(p_mat[:, ciclo, :], pm, axis=0)    
    vrv_p95 = np.percentile(v_mat[:, ciclo, :], pu, axis=0)
    prv_p95 = np.percentile(p_mat[:, ciclo, :], pu, axis=0)
    plt.plot(vrv_p5,  prv_p5,  color=cor, linestyle=':',  label=f'P{pl}')
    plt.plot(vrv_p50, prv_p50, color=cor, linestyle='-',  label=f'P{pm}')
    plt.plot(vrv_p95, prv_p95, color=cor, linestyle='--', label=f'P{pu}')
    plt.xlabel(f'{label} volume [mL]')
    plt.ylabel(f'{label} pressure [mmHg]')
    plt.legend(loc='best', bbox_to_anchor=(1, 1),frameon=False)
    plt.title(f"t={ciclo} days")
    plt.tight_layout()
    plt.savefig(f'fig_pvloop_{label}_{ciclo}.pdf',dpi=300)
    plt.savefig(f'fig_pvloop_{label}_{ciclo}.png',dpi=300)
    plt.show()    

# -----------------------------------------------------------------------------
# Seguindo o código da Camila:

def plot_cv(tf, dcv, nsa, basedir=''):

    # Sistema Cardiovascular

    # Número de arquivos em cada diretório (nciclos)
    nc = int(tf / dcv) + 1 
    # Numero de passos de tempo (nt)
    nt = 1000

    # Calcular quantidades médias para plotagem dos gráficos
    prv_mat = np.zeros((nsa,nc,nt))
    vrv_mat = np.zeros((nsa,nc,nt))
    plv_mat = np.zeros((nsa,nc,nt))
    vlv_mat = np.zeros((nsa,nc,nt))
    ppa_mat = np.zeros((nsa,nc,nt))
    vpa_mat = np.zeros((nsa,nc,nt))
    psa_mat = np.zeros((nsa,nc,nt))
    vsa_mat = np.zeros((nsa,nc,nt))
    pra_mat = np.zeros((nsa,nc,nt))
    vra_mat = np.zeros((nsa,nc,nt))

    #
    # carrega dados
    #

    # numero de amostras
    for i in range(nsa):
        nome_diretorio = basedir + f'outputs_{i:04d}'
        diretorio = os.path.join('.', nome_diretorio)
        print("processando amostra %d - %s" % (i,diretorio))
        # numero de ciclos           
        for j in range(nc):
            # nome_arquivo = f'outputs_cv_{j:02d}.txt'
            nome_arquivo = f'outputs_cv_{j:02d}.npz'
            caminho_arquivo = os.path.join(diretorio, nome_arquivo)
            #print(caminho_arquivo)

            with open(caminho_arquivo, 'r') as arquivo:
                # data = np.loadtxt(caminho_arquivo)
                zdata = np.load(caminho_arquivo)
                data = zdata['data']
                prv_mat[i,j,:] = data[:,15]
                plv_mat[i,j,:] = data[:,16]
                ppa_mat[i,j,:] = data[:,11]
                psa_mat[i,j,:] = data[:, 9]
                pra_mat[i,j,:] = data[:,13]
                vrv_mat[i,j,:] = data[:, 7]
                vlv_mat[i,j,:] = data[:, 8]
                vpa_mat[i,j,:] = data[:, 3]
                vsa_mat[i,j,:] = data[:, 1]
                vra_mat[i,j,:] = data[:, 5]

    #
    # plot dos dados - confere o amostra 0 (i)
    #
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,3))

    i = 0
    cont = 0
    days = [0,5,10,15,20,40,80,100]
    for d in days:
        j = int(d/dcv)
        s = '-'
        if(cont >= 5):
            s = '--'
        ax1.plot(vrv_mat[i,j,:], prv_mat[i,j,:], lw=1.5, label='t='+str(d), linestyle=s)
        cont = cont + 1
    j = nc-1
    ax1.set_xlabel('RV volume [mL]')
    ax1.set_ylabel('RV pressure [mmHg]')    
    # ax1.legend(loc='best', bbox_to_anchor=(1, 1),frameon=False)
    # plt.savefig('fig_pvloop_lv_0.pdf', dpi=300)
    # plt.savefig('fig_pvloop_lv_0.png', dpi=300)
    # plt.tight_layout()
    # plt.show()

    i = 0
    cont = 0

    for d in days:
        j = int(d/dcv)
        s = '-'
        if(cont >= 5):
            s = '--'
        ax2.plot(vlv_mat[i,j,:], plv_mat[i,j,:], lw=1.5, label='t='+str(d), linestyle=s)
        cont = cont + 1
    j = nc-1
    ax2.set_xlabel('LV volume [mL]')
    ax2.set_ylabel('LV pressure [mmHg]')    
    ax2.legend(loc='best', bbox_to_anchor=(1, 1),frameon=False)

    plt.savefig('fig_pvloop_rv_lv.pdf', dpi=300)
    plt.savefig('fig_pvloop_rv_lv.png', dpi=300)
    plt.tight_layout()
    plt.show()

    # sys.exit(0)

    # plot pv-loops
    plot_pvloop(vlv_mat, plv_mat, 'LV', 2, cor='r')
    plot_pvloop(vlv_mat, plv_mat, 'LV', 10, cor='r')
    plot_pvloop(vlv_mat, plv_mat, 'LV', 20, cor='r')
    
    plot_pvloop(vrv_mat, prv_mat, 'RV', 2, cor='b')
    plot_pvloop(vrv_mat, prv_mat, 'RV', 10, cor='b')
    plot_pvloop(vrv_mat, prv_mat, 'RV', 20, cor='b')

    # plot fracao de ejecao - LV
    efs_lv = np.zeros((nsa,nc))
    efs_rv = np.zeros((nsa,nc))
    for i in range(nsa):
        qoi = CalcQoIs(vlv_mat,plv_mat,vrv_mat,prv_mat,psa_mat,ppa_mat,pra_mat,i,nc)
        efs_lv[i,:] = qoi['LV_EF'][:]
        efs_rv[i,:] = qoi['RV_EF'][:]
    
    nn = int(tf/dcv) + 1
    tt = np.linspace(0,tf,nn)
    p5 = np.percentile(efs_lv, 5, axis=0)
    p95 = np.percentile(efs_lv, 95, axis=0)
    plt.plot(tt,np.mean(efs_lv, axis=0), 'k', label='mean')
    plt.fill_between(tt, p5, p95, color='0.8', alpha=0.4, label='PI90')
    #plt.plot(tt,np.percentile(efs, 5, axis=0), 'k:', label='p5')
    #plt.plot(tt,np.percentile(efs, 95, axis=0), 'k--', label='p95')
    plt.ylim([10,60])
    plt.legend(loc='lower right', frameon=False)
    plt.xlabel('time [days]')
    plt.ylabel('LV ejection fraction [%]')
    plt.title('case 2')
    plt.tight_layout()
    plt.savefig('fig_ejection_fraction_lv.pdf',dpi=300)
    plt.savefig('fig_ejection_fraction_lv.png',dpi=300)
    plt.show()

    # plot fracao de ejecao - RV
    nn = int(tf/dcv) + 1
    tt = np.linspace(0,tf,nn)
    p5 = np.percentile(efs_rv, 5, axis=0)
    p95 = np.percentile(efs_rv, 95, axis=0)
    plt.plot(tt,np.mean(efs_rv, axis=0), 'k', label='mean')
    plt.fill_between(tt, p5, p95, color='0.8', alpha=0.4, label='PI90')
    #plt.plot(tt,np.percentile(efs, 5, axis=0), 'k:', label='p5')
    #plt.plot(tt,np.percentile(efs, 95, axis=0), 'k--', label='p95')
    plt.legend(loc='lower right', frameon=False)
    plt.xlabel('time [days]')
    plt.ylabel('RV ejection fraction [%]')
    plt.title('case 2')
    plt.tight_layout()
    plt.savefig('fig_ejection_fraction_rv.pdf',dpi=300)
    plt.savefig('fig_ejection_fraction_rv.png',dpi=300)
    plt.show()

# -----------------------------------------------------------------------------    

#Função principal
if __name__ == "__main__":

    tf = 120.0  # dias
    dcv = 2.0   # dias
    nsa = 1000  # amostras
    
    # bd = 'outputs_caso1_dcv2.0_ok/'
    bd = 'outputs_caso2_dcv2.0/'

    print('basedir =',bd)
    print('samples =',nsa)
    print('tempo =',tf)
    print('deltacv =', dcv)

    plot_cv(tf, dcv, nsa, basedir=bd)
    
# -------------------------------------------------------------------------

