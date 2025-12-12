#Importando bibliotecas
import numpy as np 
import math
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('nature')

def plot_imune(arq1, out):

    # tratamento dos dados
    #dados_im = np.loadtxt(arq)
    
    # leitura npz
    zdata = np.load(arq1)
    print(zdata)
    dados_im = zdata['data']
    # dados_im = zdata['arr_0']

    print(dados_im)

    t = dados_im[0,:]
    V = dados_im[1,:]
    C = dados_im[15,:]
    IgG = dados_im[14,:]
    IgM = dados_im[13,:]

    # plots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(8, 2))

    ax1.plot(t,np.log10(V+1), lw=1.5, color='red', label='V')
    ax1.set_xlabel('time [days]')
    ax1.set_ylabel('V $\\log_{10}$[copies/mL]') # viremia
    #ax1.legend()
    # ax1.grid()

    ax2.plot(t, C, lw=1.5, color='blue', label='C')
    ax2.set_xlabel('time [days]')
    ax2.set_ylabel('C [pg/mL]') # cytokines
    #ax2.legend()
    # ax2.grid()

    ax3.plot(t, np.log2(IgG+1), lw=1.5, color='orange', label='IgG')
    ax3.set_xlabel('time [days]')
    ax3.set_ylabel('IgG $\\log_{2}$[S/CO]')
    #ax3.legend()
    # ax3.grid()

    ax4.plot(t, np.log2(IgM+1), lw=1.5, color='black', label='IgM')
    ax4.set_xlabel('time [days]')
    ax4.set_ylabel('IgM $\\log_{2}$[S/CO]')
    #ax4.legend()
    # ax4.grid()

    plt.tight_layout()
    plt.savefig(out + '.png', dpi=300)
    plt.savefig(out + '.pdf', dpi=300)
    plt.show()

def plot_imune_two_cases(arq1, arq2, out):

    # tratamento dos dados
    #dados_im = np.loadtxt(arq)
    
    # leitura npz
    zdata = np.load(arq1)
    print(zdata)
    dados_im = zdata['data']
    # dados_im = zdata['arr_0']

    print(dados_im)

    t = dados_im[0,:]
    V = dados_im[1,:]
    C = dados_im[15,:]
    IgG = dados_im[14,:]
    IgM = dados_im[13,:]

    # plots
    fig, axs = plt.subplots(2, 4, figsize=(8, 4))

    ax1, ax2, ax3, ax4 = axs[0]
    ax5, ax6, ax7, ax8 = axs[1]

    ax1.plot(t,np.log10(V+1), lw=1.5, color='red', label='V')
    ax1.set_xlabel('time [days]')
    ax1.set_ylabel('V $\\log_{10}$[copies/mL]') # viremia
    #ax1.legend()
    # ax1.grid()

    ax2.plot(t, C, lw=1.5, color='blue', label='C')
    ax2.set_xlabel('time [days]')
    ax2.set_ylabel('C [pg/mL]') # cytokines
    #ax2.legend()
    # ax2.grid()

    ax3.plot(t, np.log2(IgG+1), lw=1.5, color='orange', label='IgG')
    ax3.set_xlabel('time [days]')
    ax3.set_ylabel('IgG $\\log_{2}$[S/CO]')
    #ax3.legend()
    # ax3.grid()

    ax4.plot(t, np.log2(IgM+1), lw=1.5, color='black', label='IgM')
    ax4.set_xlabel('time [days]')
    ax4.set_ylabel('IgM $\\log_{2}$[S/CO]')
    #ax4.legend()
    # ax4.grid()

    zdata2 = np.load(arq2)
    dados_im2 = zdata2['data']

    t = dados_im2[0,:]
    V = dados_im2[1,:]
    C = dados_im2[15,:]
    IgG = dados_im2[14,:]
    IgM = dados_im2[13,:]

    ax5.plot(t,np.log10(V+1), lw=1.5, color='red', label='V')
    ax5.set_xlabel('time [days]')
    ax5.set_ylabel('V $\\log_{10}$[copies/mL]') # viremia

    ax6.plot(t, C, lw=1.5, color='blue', label='C')
    ax6.set_xlabel('time [days]')
    ax6.set_ylabel('C [pg/mL]') # cytokines

    ax7.plot(t, np.log2(IgG+1), lw=1.5, color='orange', label='IgG')
    ax7.set_xlabel('time [days]')
    ax7.set_ylabel('IgG $\\log_{2}$[S/CO]')

    ax8.plot(t, np.log2(IgM+1), lw=1.5, color='black', label='IgM')
    ax8.set_xlabel('time [days]')
    ax8.set_ylabel('IgM $\\log_{2}$[S/CO]')

    plt.tight_layout()
    plt.savefig(out + '.png', dpi=300)
    plt.savefig(out + '.pdf', dpi=300)
    plt.show()    

def plot_imune_c_cmax(arq1, arq2, out):

    Cmax = 12.0 # pg/mL
    
    # leitura npz
    zdata = np.load(arq1)
    print(zdata)
    dados_im = zdata['data']
    # dados_im = zdata['arr_0']

    print(dados_im)

    t = dados_im[0,:]
    C = dados_im[15,:]
   
    # plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 2))

    ax1.plot(t,C/Cmax, lw=1.5, color='blue', label='V')
    ax1.set_xlabel('time [days]')
    ax1.set_ylabel(r'$C/C_{max}$ [-]') # viremia

    zdata2 = np.load(arq2)
    dados_im2 = zdata2['data']

    t = dados_im2[0,:]
    V = dados_im2[1,:]
    C = dados_im2[15,:]

    ax2.plot(t, C/Cmax, lw=1.5, color='blue', label='C')
    ax2.set_xlabel('time [days]')
    ax2.set_ylabel(r'$C/C_{max}$ [-]') # cytokines
  
    plt.tight_layout()
    plt.savefig(out + '.png', dpi=300)
    plt.savefig(out + '.pdf', dpi=300)
    plt.show()        

if __name__ == "__main__":

    # arq = "outputs_0000_caso1_survivors/outputs_imune.npz"
    # out = 'fig_imune_acoplamento_caso1'
    # plot_imune(arq, out)

    # arq = "outputs_0000_caso2_non-survivors/outputs_imune.npz"
    # out = 'fig_imune_acoplamento_caso2'
    # plot_imune(arq, out)

    # arq1 = "outputs_0000_caso1_survivors/outputs_imune.npz"
    # arq2 = "outputs_0000_caso2_non-survivors/outputs_imune.npz"
    # out = 'fig_imune_acoplamento_casos'
    # plot_imune_two_cases(arq1, arq2, out)

    arq1 = "outputs_0000_caso1_survivors/outputs_imune.npz"
    arq2 = "outputs_0000_caso2_non-survivors/outputs_imune.npz"
    out = 'fig_imune_c_cmax'
    plot_imune_c_cmax(arq1, arq2, out)