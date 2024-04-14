# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 19:46:46 2023

@author: ULRICH_LUMENDO
"""
import numpy as np
from numpy import zeros, exp, pi,  ones, r_, conj, array , angle, diag, asmatrix, asarray, nonzero, arange, flatnonzero as find
from scipy.sparse import issparse, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
import importlib
import pandas as pad
import time
from variable import *
import matplotlib.pyplot as plt
import math                                                


def ordre(bus, gen, branch, n_bus, n_gen, n_br):
    """C'est pour ordonner les Bus_type pour que les calculs soient simplifiée"""
    
    brut = bus[:, bus_t].astype(int)   
    analysable = zeros(max(brut) + 1 )
    analysable[brut] = arange(bus.shape[0]) + ones(n_bus).astype(int)
    g, f, t = gen[:, Gbus].astype(int), branch[:, Fbus].astype(int), branch[:, Tbus].astype(int)
    bus[:, bus_t]    = analysable[ bus[:, bus_t].astype(int) ]
    gen[:, Gbus], branch[:, Fbus], branch[:, Tbus]  = analysable[ g],analysable[ f], analysable[ t]                     # Tableau prêt pour les calculs
    g, f, t = analysable[ g],analysable[ f], analysable[ t]
    return bus, gen, branch, (g-ones(n_gen)).astype(int), (f-ones(n_br)).astype(int), (t -ones(n_br)).astype(int)                                                   

def Spuissance(baseMVA, bus, gen, g, t, f, n_bus, n_gen, n_br):
    
    gen_busP, gen_busQ = zeros(n_bus), zeros(n_bus)        
    gen_busP[g], gen_busQ[g] = gen[:, Pgen], gen[:, Qgen]
    Sbus = ( gen_busP + 1j * gen_busQ - bus[:, Pl] - 1j * bus[:, Ql] ) / baseMVA
    
    Pd, Qd, Pg, Qg = ones(n_bus), ones(n_bus), ones(n_bus), ones(n_bus)
    i, j, k, l = nonzero(bus[:, Pl]), nonzero(bus[:, Ql]), nonzero(gen[:, Pgen]), nonzero(gen[:, Qgen])
    Pd[i], Qd[j], Pg[k], Qg[l] = zeros(n_bus)[i], zeros(n_bus)[j], zeros(n_bus)[k], zeros(n_bus)[l]
    pd, qd, pg, qg = bus[:, Pl], bus[:, Ql], gen_busP, gen_busQ
    return Sbus, Pd, Qd, Pg, Qg, pd, qd, pg, qg

def dSpuissance_dV(Ybus, V):
    I = range(len(V))
    Ibus = Ybus * V
    diagV = csr_matrix((V, (I, I)))
    diagI = csr_matrix((Ibus, (I, I)))
    diagVunit = csr_matrix((V / abs(V), (I, I)))
    dSdVm = diagV * conj(Ybus * diagVunit) + conj(diagI) * diagVunit

    dSdVa = 1j * diagV * conj(diagI - Ybus * diagV) 

    return dSdVm, dSdVa

def admitt(baseMVA, bus, branch, gen_bus , f_bus, t_bus, n_bus, n_gen, n_br):                  
    Yl = 1 / (branch[:, R] + 1j * branch[:, X])  
    Bc = branch[:, B]              
    ratio = ones(n_br)                           
    i = nonzero(branch[:, a])             
    ratio[i] = branch[i, a]               
    ratio = ratio * exp(1j * pi / 180 * branch[:, theta]) 
    
    
    
    Y11 = Yl + 1j * Bc / 2
    Y00 = Y11 / (ratio * conj(ratio))
    Y01 = - Yl / conj(ratio)
    Y10 = - Yl / ratio
    Ysh = (bus[:, Gs] + 1j * bus[:, Bs]) / baseMVA
    
    
    sup, inf, diag = {}, {}, {}
    tuplz = zip(list(f_bus), list(t_bus))
    tuples = list(tuplz)
    for j in range(max(n_bus, n_br)):
        if j < n_bus:
            diagidxf, diagidxt = list(np.where(array(f_bus) == j)[0]), list(np.where(array(t_bus) == j)[0])
            diag[(j, j)] = sum(Y00[diagidxf]) + sum(Y11[diagidxt])
        if j < n_br and tuples[j] in sup.keys():
            sup[tuples[j]] += Y01[j]
            inf[tuples[j][::-1]] += Y10[j]
        elif j < n_br and tuples[j] not in sup.keys():
            sup[tuples[j]], inf[tuples[j][::-1]] = Y01[j], Y10[j]
    
    idx = list(sup) + list(inf) + list(diag)
    lign, col = zip(*idx)
    Valeur = list(sup.values()) + list(inf.values()) + list(diag.values())
    Valeur, lign, col = array(Valeur), array(lign), array(col)
    Ybus = csr_matrix((Valeur, (lign, col)), (n_bus, n_bus))  + csr_matrix((Ysh, (range(n_bus), range(n_bus))), (n_bus, n_bus))

    return Ybus


def methode(Ybus, Sbus, x0, ref, PV, PQ):
    i = 0
    V = x0
    Va, Vm = angle(V), abs(V)
    pvpq = r_[PV, PQ]
    npv, npq = len(PV), len(PQ)

    fcomplex = V * conj(Ybus * V) - Sbus
    f = r_[  fcomplex[PV].real,
             fcomplex[PQ].real,
             fcomplex[PQ].imag  ]
    max_it  = 10**-5
    itera = 10
    
    while (True):

        if i == itera:
            return [], i
        i = i + 1
        dS_dVm, dS_dVa = dSpuissance_dV(Ybus, V)
        J00 = dS_dVa[array([pvpq]).T, pvpq].real
        J01 = dS_dVm[array([pvpq]).T, PQ].real
        J10 = dS_dVa[array([PQ]).T, pvpq].imag
        J11 = dS_dVm[array([PQ]).T, PQ].imag
        J = vstack([hstack([J00, J01]), hstack([J10, J11])], format="csr")
        Jf = -1 * spsolve(J, f)
        Vancien = V
        Va[PV] = Va[PV] + Jf[0:npv]
        Va[PQ] = Va[PQ] + Jf[npv:npv + npq]
        Vm[PQ] = Vm[PQ] + Jf[npv + npq:npv + 2*npq]
        V = Vm * exp(1j * Va)
        erreur = max(abs(Vancien - V))
        fcomplex = V * conj(Ybus * V) - Sbus
        f = r_[  fcomplex[PV].real,
                 fcomplex[PQ].real,
                 fcomplex[PQ].imag  ]
        if erreur < max_it or max(f) < max_it:
            break        
    return V, i
def execfp(casetype):
    case = importlib.import_module(casetype)
    donnee = case.case()
    start = time.time()
    donnee_bus = donnee["bus"]
    donnee_gen = donnee["gen"]
    donnee_br = donnee["branch"]
    n_bus, n_br, n_gen = len(donnee_bus), len(donnee_br), len(donnee_gen)
    Sbase = donnee["baseMVA"]
    
    
    donnee_bus, donnee_gen, donne_br, gen_bus, f_bus, t_bus   = ordre(donnee_bus, donnee_gen, donnee_br, n_bus, n_gen, n_br)
    Sbus, Pd, Qd, Pg, Qg, pd, qd, pg, qg = Spuissance(Sbase, donnee_bus, donnee_gen, gen_bus, f_bus, t_bus, n_bus, n_gen, n_br)
    Ybus = admitt(Sbase, donnee_bus, donnee_br, gen_bus , f_bus, t_bus, n_bus, n_gen, n_br)

    PV = find((donnee_bus[:, 1] == 2))
    PQ = find((donnee_bus[:, 1] == 1))
    REF = find((donnee_bus[:, 1] == 3))
    gbus = (donnee_gen[:, Gbus].astype(int) - ones(n_gen)).astype(int)            
    gbus = np.array(gbus)
    module = donnee_bus[:, Vm]
    argument = donnee_bus[:, Varg]


    x0  = module * exp(1j * pi/180 * argument)
    vcb = ones(x0.shape)  
    vcb[PQ] = 0     
    k = find(vcb[gbus])  
    x0[gbus[k]] = donnee_gen[k, Vg] / abs(x0[gbus[k]]) * x0[gbus[k]]
    s = time.time()
    V, itera = methode(Ybus, Sbus, x0, REF, PV, PQ) 
    end = time.time()
    if len(V) == 0:
        print('La méthode a divergé, choisissez bien les tensions initiaux')
    else:
        fig, ax = plt.subplots(figsize=(16, 8)) 
        t_vals = np.linspace(0, len(abs(V)), int(len(abs(V))))
        ax.scatter(t_vals, list(abs(V)), s=20, facecolor='C0', edgecolor='k', label=str(abs(V)))
        ax.plot(t_vals, [1.1]*len(abs(V)))
        ax.plot(t_vals, [0.9]*len(abs(V)))
        ax.set_xlabel('Noeuds i')
        ax.set_ylabel('Tension en pu')
        ax.set_title(f" - Surcharge admissIle de la tension réseau à {n_bus} noeuds")
        ax.set_ylim(0.7, 1.3)
        ax.grid()
        surch= 0
        surchn= []
        surchv= []
        for i in V:
            if abs(i) > 1.1 or abs(i) < 0.9 :
                ind = list(V).index(i)
                surch += 1
                surchn.append(ind+1)
                surchv.append(V[ind])
        print(f"Nous avons {surch} noeuds qui fonctionne au délà des limites établies")
        I = Ybus * V
        S = V * I.conjugate()
        Sl = S.real
        Sg = S.imag
        pgen = (Sl*Pg*Sbase) + pg + pd*Pg
        pchar = abs(Sl*Pd*Sbase) + pd - pgen*Qd
        ip = nonzero(pchar.round(1))
        pchar0 = zeros(n_bus)
        pchar0[ip] = ones(n_bus)[ip]
    
        ip = nonzero(pgen.round(1))
        pgen0 = zeros(n_bus)
        pgen0[ip] = ones(n_bus)[ip]
        
        notes = {
                 "Vm": pad.Series(abs(V).round(3), index=list(range(1, n_bus+1))),
                 "Va": pad.Series((angle(V)*180/math.pi).round(2), index=list(range(1, n_bus+1))),
                 "P gen": pad.Series(pgen.round(3), index=list(range(1, n_bus+1))),
                 "Q gen" : pad.Series(((Sg*Qg*Sbase + qg + qd*Qg)*pgen0).round(3), index=list(range(1, n_bus+1))),
                 "P char": pad.Series(pchar.round(3), index=list(range(1, n_bus+1))),
                 "Q char" : pad.Series(((abs(Sg*Qd*Sbase) + qd + qg*Qd)*pchar0).round(3), index=list(range(1, n_bus+1))),
            }
        surchf = {
                 "Noe en surch": pad.Series(surchn, index=list(range(1, len(surchn)+1))),
                 "Valeurs": pad.Series(abs(np.array(surchv)).round(2), index=list(range(1, len(surchn)+1))),
                 
            }
        
        df = pad.DataFrame(notes)
        print(df)
        if surch != 0:   
            sc = pad.DataFrame(surchf)
            print(sc)
        #print('Avec ', perte, 'comme perte de ligne du réseau')
        print("calculé en ",\
              round(end - start, 3), "sécondes avec ", itera - 1,\
                 "itérations en newton raphson")