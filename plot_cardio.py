#Importando bibliotecas
import numpy as np 
import math
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('nature')

def plot_psa_pressures(tf, dcv, casepath, cycle=0):
    
    # PV - loop (left)
    nc = int(tf / dcv)
    print(nc)
    c = cycle #0 #nc -1

    zdata = np.load(casepath + 'outputs_cv_%02d.npz' % c)
    data = zdata['data']

    t = data[:,0]
    psa  = data[:,9]

    plt.plot(t, psa, 'k', lw=2)
    plt.xlabel('time [mL]')
    plt.ylabel(r'$P_{sa}$ [mmHg]')
    plt.xlim([20,30])
    plt.tight_layout()
    # plt.legend()
    plt.savefig(casepath + f'fig_cv_psa_{c}.pdf',dpi=300)
    plt.show()

def plot_pressures(tf, dcv, casepath, cycle=0):
    
    # PV - loop (left)
    nc = int(tf / dcv)
    print(nc)
    c = cycle #0 #nc -1

    zdata = np.load(casepath + 'outputs_cv_%02d.npz' % c)
    data = zdata['data']

    t = data[:,0]
    prv  = data[:,15]
    plv  = data[:,16]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6,2))
    #fig.suptitle('Horizontally stacked subplots')
    ax1.plot(t, prv, 'b', lw=1.5)
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel(r'$P_{rv}$ [mmHg]')
    ax1.set_xlim([20,30])

    ax2.plot(t, plv, 'r', lw=1.5)
    ax2.set_xlabel('time [s]')
    ax2.set_ylabel(r'$P_{lv}$ [mmHg]')
    ax2.set_xlim([20,30])
    ax2.set_ylim([0,175])
    # plt.plot(t, plv, 'k', lw=2)

    plt.tight_layout()
    # plt.legend()
    plt.savefig(casepath + f'fig_cv_prv_plv_{c}.pdf',dpi=300)
    plt.show()

def plot_prv_pressures(tf, dcv, casepath, cycle=0):
    
    # PV - loop (left)
    nc = int(tf / dcv)
    print(nc)
    c = cycle #0 #nc -1

    zdata = np.load(casepath + 'outputs_cv_%02d.npz' % c)
    data = zdata['data']

    t = data[:,0]
    psa  = data[:,15]

    plt.plot(t, psa, 'k', lw=1)
    plt.xlabel('time [mL]')
    plt.ylabel(r'$P_{rv}$ [mmHg]')
    plt.xlim([20,30])
    
    plt.tight_layout()
    # plt.legend()
    plt.savefig(casepath + f'fig_cv_prv_{c}.pdf',dpi=300)
    plt.show()

def plot_pv_loop(tf, dcv, casepath):
    
    nc = int(tf / dcv)
    for i in range(nc+1):
        zdata = np.load(casepath + 'outputs_cv_%02d.npz' % i)
        data = zdata['data']

        Vlv = data[:,8]
        Plv = data[:,16]

        if(i==0):
            plt.plot(Vlv, Plv, label ='Paciente normal') 
        else:
            plt.plot(Vlv, Plv, label ='Janela %02d' % (i))

    plt.xlabel('Vlv [mL]')
    plt.ylabel('Plv [mmHg]')
    plt.tight_layout()
    plt.legend()
    plt.title('Pressão x Volume do ventrículo esquerdo')
    plt.savefig(casepath + 'fig_pv_loop_left.pdf',dpi=300)
    plt.show()

if __name__ == "__main__":

    # -------------------------------------------------------------------------

    tf = 120.0 #dias
    dcv = 2.0 #dias

    # -------------------------------------------------------------------------

    casepath = 'outputs_fig4/'

    # plot_pv_loop(tf, dcv, casepath)
    
    # plot_psa_pressures(tf, dcv, casepath, cycle=0)
    
    plot_pressures(tf, dcv, casepath, cycle=0)
