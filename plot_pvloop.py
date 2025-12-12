#Importando bibliotecas
import numpy as np 
import math
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('nature')

def plot_pv_loop(tf, dcv, cycles, casepath):
    
    nc = int(tf / dcv)
    # for i in range(0,nc+1,6):
    for i in cycles:
        zdata = np.load(casepath + 'outputs_cv_%02d.npz' % i)
        data = zdata['data']
        Vlv = data[:,8]
        Plv = data[:,16]
        print(i)
        if(i==0):
            plt.plot(Vlv, Plv, lw=1.5, label ='initial') 
        elif(i<13):
            plt.plot(Vlv, Plv, lw=1.5, label ='%02d' % (i)) 
        else:
            plt.plot(Vlv, Plv, '--', label ='%02d' % (i))

    plt.xlabel('Vlv [mL]')
    plt.ylabel('Plv [mmHg]')
    # plt.legend()
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.tight_layout()
    plt.savefig(casepath + 'fig_pvloop_left.pdf',dpi=300)
    plt.show()

if __name__ == "__main__":

    # -------------------------------------------------------------------------

    tf = 120.0 #dias
    dcv = 2.0 #dias

    # -------------------------------------------------------------------------

    casepath = 'outputs_caso1_dcv2.0/outputs_0000/'
    ciclos = [0,1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50,60]

    plot_pv_loop(tf, dcv, ciclos, casepath)
    