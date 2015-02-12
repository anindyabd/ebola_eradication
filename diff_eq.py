__author__ = "aguha@colgate.edu" 

from numpy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages

def deriv(vector, t, beta_I, beta_H, beta_F, alpha, gamma_H, gamma_I, gamma_D, gamma_DH, gamma_F, gamma_IH, delta1, delta2, delta3, iota):
    
    S = vector[0]
    E = vector[1]
    I = vector[2]
    H = vector[3]
    F = vector[4]
    M = vector[5]
    R = vector[6]

    N = S + E + I + H + F + M + R 
    
    Sprime =  - ( beta_I * S * I + beta_H * S * H + beta_F * S * F ) / N 
    Eprime = ( ( beta_I * S * I + beta_H * S * H + beta_F * S * F ) / N ) - alpha * E
    Iprime = alpha * E - ( gamma_H * iota + gamma_I * (1 - iota) * (1 - delta1) + gamma_D * (1 - iota) * delta1 ) * I 
    Hprime = gamma_H * iota * I - ( gamma_DH * (delta2 + delta3)  + gamma_IH * (1 - (delta2 + delta3)) ) * H 
    Fprime = gamma_D * (1 - iota) * delta1 * I + gamma_DH * (delta2 + delta3) * H - gamma_F * F 
    Mprime = gamma_IH*(1 - (delta2 + delta3))*H  
    Rprime = gamma_I * (1 - iota) * (1 - delta1) * I + gamma_F*F 

    return [Sprime, Eprime, Iprime, Hprime, Fprime, Mprime, Rprime]

def run_odeint(beta_I=0.16, beta_H=0.062, beta_F=0.489, alpha=(1/12.0), gamma_H=(1/3.24), gamma_I=(1/15.0), gamma_D=(1/13.31), gamma_DH=(1/10.07), gamma_F=(1/2.01), gamma_IH=(1/15.88), delta1=0.5, delta2=0.5, delta3=0, iota=0.197): 
    
    init = [1667843, 20, 131, 162, 0, 142, 2000]
    times = linspace(0, 300, 400)
    
    y = odeint(deriv, init, times, args=(beta_I, beta_H, beta_F, alpha, gamma_H, gamma_I, gamma_D, gamma_DH, gamma_F, gamma_IH, delta1, delta2, delta3, iota))

    return y 

def draw_figure(y, figno, figname):

    times = linspace(0, 300, 400)
    plt.figure(figno)
    plt.plot(times, y[:,0], '-g', label='S') 
    plt.plot(times, y[:,1], '-y', label='E')
    plt.plot(times, y[:,2], '-b', label='I') 
    plt.plot(times, y[:,3], '-m', label='H') 
    plt.plot(times, y[:,4], '-c', label='F')
    plt.plot(times, y[:,5], '#2E0854', label='M')  
    plt.plot(times, y[:,6], '-r', label='R') 
    plt.legend() 
    #plt.title('Figure ' + str(figno))
    plt.xlabel('Time (days)')
    plt.ylabel('Population')
    plt.savefig(figname)
    #plt.hold(False)

if __name__ == '__main__':


    draw_figure(run_odeint(), 1, 'no_intervention.png')
    draw_figure(run_odeint(delta2=0.2, delta3=0.25), 2, 'vaccination_early_stage.png')
    draw_figure(run_odeint(delta2=0.1, delta3=0.25), 3, 'potent_vaccination_early_stage.png')
    draw_figure(run_odeint(delta2=0.05, delta3=0.25), 4, 'more_potent_vaccination_early_stage.png')
    draw_figure(run_odeint(delta2=0, delta3=0.25), 5, 'fully_potent_vaccination_early_stage.png')
    
    draw_figure(run_odeint(delta2=0.1, delta3=0.1), 6, 'fully_potent_vaccination_all_stages.png')

    draw_figure(run_odeint(beta_H=0.0375, delta2=0.0, delta3=0.25, iota=0.24625), 7, 'identify_and_isolate_es_1.png')
    draw_figure(run_odeint(beta_H=0.025, delta2=0.0, delta3=0.25, iota=0.2955), 8, 'identify_and_isolate2_es_2.png')
    draw_figure(run_odeint(beta_H=0.00125, delta2=0.0, delta3=0.25, iota=0.34475), 9, 'identify_and_isolate3_es_3.png')
   
    draw_figure(run_odeint(beta_I=0.12, beta_H=0.0375, delta2=0, delta3=0.25, iota=0.24625, gamma_F=(1/1.51)), 10, 'close_down1_es.png')
    draw_figure(run_odeint(beta_I=0.11, beta_H=0.025, delta2=0, delta3=0.25, iota=0.2955, gamma_F=(1/1.21)), 11, 'close_down2_es.png')
    
