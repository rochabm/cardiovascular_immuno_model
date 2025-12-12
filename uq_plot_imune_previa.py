import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import*
import pylab as pl
import os
import sys


def plot_imune():

    #nsa -> nº de samples
    nsa = 15  # Número de amostras
    # -------------------------------------------------------------------------
    tf = 120.0 #dias
    dcv = 5.0 #dias
    nc = int(tf / dcv) + 1
    # -------------------------------------------------------------------------
    #Com nt também?
    V_mat = np.zeros((nsa,nc))
    C_mat = np.zeros((nsa,nc))
    IgG_mat = np.zeros((nsa,nc))
    IgM_mat = np.zeros((nsa,nc))

    # numero de amostras
    for i in range(nsa):
        nome_diretorio = f'outputs_imune_{i:04d}'
        diretorio = os.path.join('.', nome_diretorio)
        nome_arquivo = f'outputs_imune.txt'
        caminho_arquivo = os.path.join(diretorio, nome_arquivo)
        with open(caminho_arquivo, 'r') as arquivo:
            dados_im = np.loadtxt(caminho_arquivo)
            t = dados_im[0,:]
            V = dados_im[1,:]
            C = dados_im[15,:]
            IgG = dados_im[14,:]
            IgM = dados_im[13,:]

            V_mat[i,:] = V
            C_mat[i,:] = C
            IgG_mat[i,:] = IgG
            IgM_mat[i,:] = IgM

    #Vetores para armazenamento de médias e desvios padrão
    V_mean = np.zeros(nc)
    V_std = np.zeros(nc)
    C_mean = np.zeros(nc)
    C_std = np.zeros(nc)
    IgG_mean = np.zeros(nc)
    IgG_std = np.zeros(nc)
    IgM_mean = np.zeros(nc)
    IgM_std = np.zeros(nc)

    #Laço que calcula e armazena para cada ciclo nc
    for j in range(nc):
        V_mean[j] = np.mean(V_mat[i,:])
        V_std[j] = np.std(V_mat[i,:])
        C_mean[j] = np.mean(C_mat[i,:])
        C_std[j]= np.std(C_mat[i,:])
        IgG_mean[j] = np.mean(IgG_mat[i,:])
        IgG_std[j] = np.std(IgG_mat[i,:])
        IgM_mean[j] = np.mean(IgM_mat[i,:])
        IgM_std[j] = np.std(IgM_mat[i,:])


    #Plotagens de gráfico (Ajeitar)
    fig, axs = plt.subplots(2, 2)

    # Plotar o gráfico Viremia
    axs[0, 0].plot(t[k:], np.log10(V_mean+1), color='red', label='V')

    #Conferir esses logs e afins
    axs[0, 0].fill_between(t[k:], np.log10(V_mean+1) - np.log10(V_std+1), np.log10(V_mean+1) + np.log10(V_std+1), alpha=0.25, color='darkgreen')
    axs[0, 0].set_title('Viremia')
    axs[0, 0].set(xlabel='Tempo [dias]', ylabel='Viremia $\\log_{10}$[cópias/mL]')

    # Plotar o gráfico Citocinas
    axs[0, 1].plot(t, C_mean, color='blue', label='C')
    axs[0, 1].fill_between(t, C_mean - C_std, C_mean + C_std, alpha=0.25, color='darkgreen')
    axs[0, 1].set_title('Citocinas')
    axs[0, 1].set(xlabel='Tempo [dias]', ylabel='Citocinas [pg/mL]')

    # Plotar o gráfico IgG
    axs[1, 0].plot(t, np.log2(IgG_mean+1), color='orange', label='IgG')
    axs[1, 0].fill_between(t, np.log2(IgG_mean+1) - np.log2(IgG_std+1), np.log2(IgG_mean+1) + np.log2(IgG_std+1), alpha=0.25, color='darkgreen')
    axs[1, 0].set_title('IgG')
    axs[1, 0].set(xlabel='Tempo [dias]', ylabel='IgG $\\log_{2}$[S/CO]')

    # Plotar o gráfico IgM
    axs[1, 1].plot(t, np.log2(IgM_mean+1), color='black', label='IgM')
    axs[1, 1].fill_between(t, np.log2(IgM_mean+1) - np.log2(IgM_std+1), np.log2(IgM_mean+1) + np.log2(IgM_std+1), alpha=0.25, color='darkgreen')
    axs[1, 1].set_title('IgM')
    axs[1, 1].set(xlabel='Tempo [dias]', ylabel='IgM $\\log_{2}$[S/CO]')

    plt.tight_layout()
    plt.savefig('imune.png',format = 'png',dpi=300)
    plt.show()
  
if __name__ == "__main__":
    plot_imune()

