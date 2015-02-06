__author__ = "aguha@colgate.edu" 

from numpy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages

def deriv(vector, t, beta_I, beta_H, beta_F, alpha, gamma_H, gamma_I, gamma_D, gamma_DH, gamma_F, gamma_IH, delta1, delta2, iota):
    
    S = vector[0]
    E = vector[1]
    I = vector[2]
    H = vector[3]
    F = vector[4]
    R = vector[5]

    N = S + E + I + H + F + R 
    
    Sprime = ( beta_I * S * I + beta_H * S * H + beta_F * S * F ) / N 
    Eprime = ( ( beta_I * S * I + beta_H * S * H + beta_F * S * F ) / N ) - alpha * E
    Iprime = alpha * E - ( gamma_H * iota + gamma_I * (1-iota) * (1 - delta1) + gamma_D * (1 - iota) * delta1 ) * I 
    Hprime = gamma_H * iota * I - ( gamma_DH * delta2 + gamma_IH * (1 - delta2) ) * H 
    Fprime = gamma_D * (1 - iota) * delta1 * I + gamma_DH * delta2 * H - gamma_F * F 
    Rprime = gamma_I * (1 - iota) * (1 - delta1) * I + gamma_IH * (1 - delta2) * H + gamma_F*F 

    return [Sprime, Eprime, Iprime, Hprime, Fprime, Rprime]

def run_odeint(init, times, beta_I=0.16, beta_H=0.062, beta_F=0.489, alpha=(1/12), gamma_H=(1/3.24), gamma_I=(1/15), gamma_D=(1/13.31), gamma_DH=(1/10.07), gamma_F=(1/2.01), gamma_IH=(1/15.88), delta1=0.5, delta2=0.5, iota=0.197): 
    
    y = odeint(deriv, init, times, args=(beta_I, beta_H, beta_F, alpha, gamma_H, gamma_I, gamma_D, gamma_DH, gamma_F, gamma_IH, delta1, delta2, iota))

    return y 

if __name__ == '__main__':
    
    init = [5000, 1000, 100, 800, 500, 200]
    times = linspace(0, 300, 400)

    y = run_odeint(init, times) 

    plt.plot(times, y[:,0], '-r', label='S') 
    plt.plot(times, y[:,1], '-g', label='E')
    plt.plot(times, y[:,2], '-b', label='I') 
    plt.plot(times, y[:,3], '-m', label='H') 
    plt.plot(times, y[:,4], '-c', label='F')  
    plt.plot(times, y[:,5], '-y', label='R') 
    plt.legend() 

    plt.show()

