#!/usr/bin/env python
#
# numerical iterations of evolving elevation, 
# a, plasticity, b, and maternal effects, m.
#
# This work is licensed under a Creative Commons 
# Attribution-NonCommercial-ShareAlike 4.0 
# International License. 
# http://creativecommons.org/licenses/by-nc-sa/4.0/
#
# when using or extending the code, please cite
# Kuijper, B & Hoyle, R. B. 2015
# When to rely on maternal effects and when to rely on phenotypic plasticity?
# Evolution, http://dx.doi.org/10.1111/evo.12635 


import os, random, math, copy, sys
import numpy as np
from datetime import date, datetime


class IterM:

    # skip in data output
    skip = 1

    # stochastic autocorr
    rho = 0.5

    # mean of stationary autocorrelated noise
    xi_bar = 0; 

    # standard deviation of stationary autocorrelated noise
    sig_xi = 0.01

    # strength of selection on z
    omega_z_2 = 40

    # variance in environmental component
    sigma_e2 = 1.0 

    # intercept of reference envt
    A=0;
    B=2;

    ampl = 1

    # maximal fitness
    Wmax = 1.0

    # initial value of the envt
    epsinit = 0

    # final value of the environment
    epsend = 10

    # length of simulation
    nrp = 500000

    Uton = 50000

    Utoff = 5000000

    # initialize a single simulation run with particular values
    # for Gaa, Gbb, Gmm, omega_b_2, omega_m_2 which override the class
    # defaults specified above
    # shift: logical variable stating whether we are interested in a single environmental
    # shift (True) or in a sinusoidal environment (False)
    # freq: frequency of the sinusoidal
    # tau: the delay in environmental information
    # 

    def __init__(self, Gaa, Gbb, Gmm, omega_b_2, omega_m_2, shift=True, freq=0, tau=0.25, delta=10, init_g=0.001,init_b=0.001,init_m=0.001):

        self.Gaa = Gaa
        self.Gbb = Gbb
        self.Gmm = Gmm
        self.omega_b_2 = omega_b_2
        self.omega_m_2 = omega_m_2
        self.freq = freq
        self.tau = tau

        # determine shape of environmental shift
        # baseline environment
        self.delta = delta

        self.init_m = init_m
        self.init_g = init_g
        self.init_b = init_b

        print("\n\nGaa:;" + str(self.Gaa) + "\n" + \
                "Gbb:;" + str(self.Gbb) + "\n" + \
                "Gmm:;" + str(self.Gmm) + "\n" + \
                "delta:;" + str(self.delta) + "\n" + \
                "B:;" + str(self.B) + "\n" + \
                "init_g:;" + str(self.init_g) + "\n" + \
                "init_b:;" + str(self.init_b) + "\n" + \
                "init_m:;" + str(self.init_m) + "\n" + \
                "omega_b_2:;" + str(self.omega_b_2) + "\n" + \
                "omega_m_2:;" + str(self.omega_m_2) + "\n" + \
                "omega_z_2:;" + str(self.omega_z_2) + "\n" + \
                "freq:;" + str(self.freq) + "\n" + \
                "tau:;" + str(self.tau))

        # storage for variables
        vec = [0 for i in range(0,self.nrp)] 

        # generate empty vectors 
        # to store the data
        self.abar = copy.deepcopy(vec)
        self.bbar = copy.deepcopy(vec)
        self.mbar = copy.deepcopy(vec)
        self.zbar = copy.deepcopy(vec)
        self.zsbar = copy.deepcopy(vec)
        self.gaz = copy.deepcopy(vec)
        self.gbz = copy.deepcopy(vec)
        self.gmzstar = copy.deepcopy(vec)
        self.fitness = copy.deepcopy(vec)
        self.epst = copy.deepcopy(vec)
        self.epstT = copy.deepcopy(vec)
        self.epstT1 = copy.deepcopy(vec)
        self.epstT2 = copy.deepcopy(vec)
        self.sigz = copy.deepcopy(vec)
        self.xi_devel = copy.deepcopy(vec)
        self.xi_select = copy.deepcopy(vec)

        # initialization
        self.mbar[0] = self.mbar[1] = self.mbar[2]  = self.init_m
        self.abar[0] = self.abar[1] = self.abar[2]  = self.init_g
        self.bbar[0] = self.bbar[1] = self.bbar[2]  = self.init_b
        self.zsbar[0] = self.zsbar[1]  = 0.0
        self.zbar[0] = self.zbar[1]  = 0.0
        self.gmzstar[0] = self.gmzstar[1] = 0.0

        if self.mbar[0] == 2:
            self.gaz[0] = self.gaz[1] = 0.0
        else:
            self.gaz[0] = self.gaz[1] = self.Gaa / (2 - self.mbar[0])


        self.gbz[0] = self.gbz[1] = 0
        self.epstT[0] = self.epstT[1] = self.epsinit
        self.epst[0] = self.epst[1] = self.epsinit
        self.GG = [ [ self.Gaa, 0, 0], [ 0, self.Gbb, 0] , [ 0, 0, self.Gmm] ]

        # set the type of environmental change for the current simulation
        self.envt_fun = self.envt_shift
        self.envtT_fun = self.envtT_shift

        if not shift:
            self.envt_fun = self.envt_sin
            self.envtT_fun = self.envtT_sin

        # store for the write out
        self.shift = shift

    # environment changes in a single shift
    # t: the current time t
    def envt_shift(self, t):
        return(self.epsinit + self.delta*(t>self.Uton)*(t<self.Utoff) + (self.epsend - self.epsinit) * (t>=self.Utoff))

    def envtT_shift(self, t):
        return(self.epsinit + self.delta*((t-self.tau)>self.Uton)*((t-self.tau)<self.Utoff) + (self.epsend - self.epsinit) * ((t-self.tau)>=self.Utoff))

    # environment changes periodically 
    # t: the current time t
    def envt_sin(self, t):
        return(self.epsinit+self.ampl*math.sin(self.freq *t))

    def envtT_sin(self, t):
        return(self.epsinit+self.ampl*math.sin(self.freq * (t-self.tau)))

    # run the actual simulation
    def run(self):

        # generate environmental sequence
        self.xi_devel[0] = self.sig_xi * random.gauss(0, 1);

        for k in range(0,self.nrp-1):
            if (self.rho == 0.0):
                self.xi_select[k] = self.sig_xi * random.gauss(0,1);
                self.xi_devel[k] = self.sig_xi * random.gauss(0,1);
            else:
                tmp = (self.rho*self.xi_devel[k]) + (self.sig_xi * math.sqrt(1-(self.rho**2))*random.gauss(0,1))
                self.xi_select[k]=tmp

                if (self.tau > 0):
                    j=1.0/self.tau

                    # get through timesteps until you get the actual developmental environment
                    if j >= 2:
                        for i in range(0,int(round(j))):
                            tmp =(self.rho*tmp) + (self.sig_xi * math.sqrt(1-(self.rho**2))*random.gauss(0,1))

                        self.xi_devel[k+1] = (self.rho*tmp) + (self.sig_xi * math.sqrt(1-(self.rho**2))*random.gauss(0,1))

        # the guts of the code: loop through time
        for t in range(2,self.nrp-2):

            # environmental change
            epsfun = self.envt_fun(t)
            eps_t = epsfun + self.xi_select[t]
            self.epst[t] = eps_t 
            
            epstTfun = self.envtT_fun(t)
            epstT = epstTfun + self.xi_devel[t]
            self.epstT[t] = epstT

            epstT1 = self.epstT[t-1] + self.xi_devel[t-1]
            epstT2 = self.epstT[t-2] + self.xi_devel[t-2]

            sigma_z2 = self.Gaa + self.Gbb*epstT**2.0 + self.sigma_e2 \
                    + 2 * self.mbar[t] * (self.gaz[t-1] + epstT * self.gbz[t-1])

            sigma_z2 += self.gmzstar[t-1]**2 \
                    + 2*self.gmzstar[t-1]*self.mbar[t]*self.zsbar[t-1] \
                    + self.Gmm * self.zsbar[t-1]**2

            sigma_z2 = sigma_z2 / (1.0 - self.Gmm - self.mbar[t]**2)

            if sigma_z2 < 0:
                print("falen: sigma_neg")

            sigma_z2 = max(sigma_z2,0)
            self.sigz[t] = sigma_z2

            gamma = 1.0 / (self.omega_z_2 + sigma_z2)
            gamma_b = 1.0 / (self.omega_b_2 + self.Gbb)
            gamma_m = 1.0 / (self.omega_m_2 + self.Gmm)

            # update \bar{z}_{t}
            self.zbar[t] = self.abar[t] + (self.bbar[t]*epstT) + self.gmzstar[t-1] + self.mbar[t]*self.zsbar[t-1]

            # update the value of theta
            theta = self.A + self.B * eps_t

            tmp = self.zbar[t] - theta

            dzbara =  1+self.mbar[t]
            dzbarb = epstT + self.mbar[t] * epstT1
            dzbarm = self.zsbar[t-1] + self.mbar[t]*self.zsbar[t-2] + 0.5*self.gmzstar[t-1]

            dsiga = 2 * (self.mbar[t]*self.gmzstar[t-1] + self.Gmm*self.zsbar[t-1])/(1.0 - self.Gmm - self.mbar[t]**2)

            dsigb = dsiga * epstT1

            dsigm = (2.0 / (1.0 - self.Gmm - self.mbar[t]**2)) * (self.gaz[t-1] * (1.0 + 0.5 * self.mbar[t]) + (epstT + 0.5* self.mbar[t] * epstT1) * self.gbz[t-1] + 0.5 * self.gmzstar[t-1]**2 + self.gmzstar[t-1] * self.zsbar[t-1] * (1.0 + 0.5 * self.mbar[t]) + self.gmzstar[t-1] * self.mbar[t] * self.zsbar[t-2] + self.Gmm * self.zsbar[t-1] * self.zsbar[t-2]) + 2 * self.mbar[t] * sigma_z2 / (1-self.Gmm-self.mbar[t]**2)

            beta = np.array([ (1.0 / self.omega_z_2) * ( - tmp * dzbara - .5 * dsiga),\
                (1.0 / self.omega_z_2) * (-tmp * dzbarb - .5 * dsigb - self.omega_z_2 * self.bbar[t] / self.omega_b_2),\
                (1.0 / self.omega_z_2) * (-tmp * dzbarm - .5 * dsigm - self.omega_z_2 * self.mbar[t] / self.omega_m_2)])

            self.GG = np.array([[self.Gaa, 0.0, 0.0],[0.0, self.Gbb, 0.0], [ 0.0, 0.0, self.Gmm ]])

            delta_change = np.dot(self.GG,beta)

            prefacfit = math.sqrt(gamma * gamma_b * gamma_m * self.omega_z_2 * self.omega_b_2 * self.omega_m_2)

            self.fitness[t] = self.Wmax * prefacfit * math.exp(-(gamma/2.0) * ((self.zbar[t] - theta)**2) - .5 * gamma_m * (self.mbar[t]**2) - .5 * gamma_b * (self.bbar[t]**2))

            self.abar[t+1] = self.abar[t] + delta_change[0]
            self.bbar[t+1] = self.bbar[t] + delta_change[1]
            self.mbar[t+1] = self.mbar[t] + delta_change[2]

            self.zsbar[t] = self.abar[t+1] + self.bbar[t+1] * epstT + self.gmzstar[t-1] + self.mbar[t+1] * self.zsbar[t-1]

            # update C_{m_{t}z_{t-1}^*}
            self.gmzstar[t] = 0.5 * self.mbar[t+1] * self.gmzstar[t-1] + 0.5 * self.Gmm * self.zsbar[t-1]

            self.gaz[t] = 0.5 * self.Gaa + 0.5 * self.mbar[t+1] * self.gaz[t-1]

            self.gbz[t] = 0.5 * self.Gbb*epstT + 0.5 * self.mbar[t+1] * self.gbz[t-1]


    def printnums(self):

        d = datetime.now()
        filename = "iter_" + d.strftime("%d_%m_%Y_%H%M%S_%f") + ".csv"

        a = open(filename,"w")

        a.write("time;abar;bbar;mbar;zbar;zbartmin1;envt;fitness;soft;\n")

        num_round = 5

        # loop over the data and write it out
        for i in range(0,len(self.mbar),1):

            soft = 0

            if i > 0 and self.zbar[i] > 0:
                soft = self.mbar[i] * self.bbar[i] * self.epst[i-1] / (self.zbar[i])

            if i < 100000 or i % self.skip == 0:
                a.write(str(i) + ";" + str(round(self.abar[i],num_round)) + ";" \
                        + str(round(self.bbar[i],num_round)) + ";" + str(round(self.mbar[i],num_round)) + ";" \
                        + str(round(self.zbar[i],num_round)) + ";" + str(round(self.zsbar[i],num_round)) + ";" \
                        + str(round(self.epst[i],num_round)) + ";" \
                        + str(round(self.fitness[i],num_round)) + ";" + str(round(soft,num_round)) + ";\n")

        a.write("\n\nGaa:;" + str(self.Gaa) + "\n" + \
                "Gbb:;" + str(self.Gbb) + "\n" + \
                "Gmm:;" + str(self.Gmm) + "\n" + \
                "delta:;" + str(self.delta) + "\n" + \
                "B:;" + str(self.B) + "\n" + \
                "A:;" + str(self.A) + "\n" + \
                "omega_b_2:;" + str(self.omega_b_2) + "\n" + \
                "omega_m_2:;" + str(self.omega_m_2) + "\n" + \
                "omega_z_2:;" + str(self.omega_z_2) + "\n" + \
                "shift:;" + str(self.shift) + "\n" + \
                "sig_xi:;" + str(self.sig_xi) + "\n" + \
                "rho:;" + str(self.rho) + "\n" + \
                "freq:;" + str(self.freq) + "\n" + \
                "tau:;" + str(self.tau)) 

        a.close()

# end class 

Gaa = float(sys.argv[1])
Gbb = float(sys.argv[2])
Gmm = float(sys.argv[3])
omega_b_2 = float(sys.argv[4])
omega_m_2 = float(sys.argv[5])
shift = bool(int(sys.argv[6]))
freq = float(sys.argv[7])
tau = float(sys.argv[8])
delta = float(sys.argv[9])
init_g = float(sys.argv[10])
init_b = float(sys.argv[11])
init_m = float(sys.argv[12])

a = IterM(Gaa, Gbb, Gmm, omega_b_2, omega_m_2, shift, freq, tau, delta, init_g, init_b, init_m)
a.run()
a.printnums()
