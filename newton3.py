# -*- coding: utf-8 -*-

"""
Created on Sat Sep 23 21:11:23 2023
@author: ULRICH_LUMENDO
"""
import time
from math import sin, cos, pi
import numpy as np
import sympy as sp
from sys import stderr
from numpy import ones, conj, nonzero, any, exp, pi, r_
from scipy.sparse import csr_matrix
import cmath
start1 = time.time()



from case import case




unknows_part = {}
donnee = case()
donnee_bus = donnee["bus"]
donnee_gen = donnee["gen"]
donnee_br = donnee["branch"]
nbr = len(donnee_bus)
base = donnee["baseMVA"]


def puissance(nbr, input_data, admittance, ref, pq, slack):
    difff = []      #liste contenant les fonctions de la puissance réelle
    diffg = []      #liste contenant les fonctions de la puissance apparente
    difffv = []      #liste contenant les fonctions de la puissance réelle
    diffgv = []      #liste contenant les fonctions de la puissance apparente
    

    for e in range(1, nbr+1):
        if e not in ref:
            
            #initialisation des fonctions
            f = [0, []]  #la puissance réelle
            g = [0, []]  #la puissance apparente
            
            
            #variable ou donnée en volt
            if 'v' + str(e) in input_data.keys():
                globals()['x' + str(e)] = input_data['v' + str(e)]
            else:
                globals()['x' + str(e)] = sp.symbols('x' + str(e))
                
            
            #variable ou donnée en phase
            if 'sigma' + str(e) in input_data.keys():
                globals()['z' + str(e)] = input_data['sigma' + str(e)]
            else:
                globals()['z' + str(e)] = sp.symbols('z' + str(e))
                
            #celui prénant en compte l'expression somme
            for n in range(1, nbr+1):
                if 'v' + str(n) in input_data.keys():
                    globals()['x' + str(n)] = input_data['v' + str(n)]
                else:
                    globals()['x' + str(n)] = sp.symbols('x' + str(n)) 
                if 'sigma' + str(n) in input_data.keys():
                    globals()['z' + str(n)] = input_data['sigma' + str(n)]
                else:
                    globals()['z' + str(n)] = sp.symbols('z' + str(n))
                #première élément de la somme
                if admittance[int(e)-1][int(n)-1] != 0:
                    
                    g1 = globals()['x' + str(e)] * globals()['x' + str(n)] * admittance[int(e)-1][int(n)-1].real * sp.sin(globals()['z' + str(e)] - globals()['z' + str(n)]) - globals()['x' + str(e)] * globals()['x' + str(n)] * admittance[int(e)-1][int(n)-1].imag * sp.cos(globals()['z' + str(e)] - globals()['z' + str(n)])
                
                    f1 = globals()['x' + str(e)] * globals()['x' + str(n)] * admittance[int(e)-1][int(n)-1].real * sp.cos(globals()['z' + str(e)] - globals()['z' + str(n)]) + globals()['x' + str(e)] * globals()['x' + str(n)] * admittance[int(e)-1][int(n)-1].imag * sp.sin(globals()['z' + str(e)] - globals()['z' + str(n)])
                    #incrementation dans les fonctions de départ
                    f[0] += f1
                    f[1] += [globals()['x' + str(n)], globals()['x' + str(e)], globals()['z' + str(e)], globals()['z' + str(n)]]
                    g[0] += g1
                    g[1] += [globals()['x' + str(n)], globals()['x' + str(e)], globals()['z' + str(e)], globals()['z' + str(n)]]
                    
                else:
                    g1 = 0
                    f1 = 0
                    f[0] += f1
                    g[0] += g1
            #ajout dans les listes de puissances
            f[1] = list(set(f[1]))
            g[1] = list(set(g[1]))
            
            difff.append(f[0])
            diffg.append(g[0])
            difffv.append(f[1])
            diffgv.append(g[1])
    diff = []       #initialisation de la matrice des fonctions de puissances fx
    #ajout des valeurs initiales
    compteur = 0
    codiff2 = 0
    while True:
        if compteur == nbr:
            break
        if compteur +1 not in ref:
            difff[codiff2]+= input_data['p' + str(compteur + 1)] 
            codiff2 += 1
            compteur += 1
        else:
            compteur += 1
    #ajout des valeurs initiales faire attention au swing bus
    compteur = 0
    codiff2 = 0
    diffg2 = []
    diffg2v = []
    while True:
        if compteur == nbr:
            break
        if compteur+1 in pq:
            diffg[codiff2] += input_data['q' + str(compteur + 1)]
            diffg2.append(diffg[codiff2])
            diffg2v.append(diffgv[codiff2])
            codiff2 += 1
            compteur += 1
        elif compteur+1 not in pq and compteur +1 in slack:
            codiff2 += 1
            compteur += 1
        elif compteur+1 not in pq:
            compteur += 1
        else:
            compteur += 1
    #matrice des fonctions de puissances fx
    diff1 = difff + diffg2
    diff1v = difffv + diffg2v
    #calcul des différentialles totales et ajout dans la matrice jordan jx
    return diff1, diff1v




def Ybus(donnee_bus, donnee_gen, donnee_br, nbr):
    input_data = {}                 #dictionnaire des donnÃ©es
    unknows = []                    #listes des inconnues
    slack = []
    pq = []
    ref = []
    #rÃ©cuperation de type de bus
    v=0
    s=0
    j = 0
    compteur = 0
    for i in donnee_bus:   
        if i[1] == 1:
            unknows.append('sigma' + str(j+1))
            unknows.append('v' + str(j+1))
            
            donneep = 0
            donneeq = 0
            index = 0
            donnee_genu = []
            for i in range(len(donnee_gen)):
                donnee_genu.append(donnee_gen[i][0])
            for valeur in donnee_genu:
                if valeur == j+1:
                    donneep += donnee_gen[index][1]
                    donneeq += donnee_gen[index][2]
                index += 1
            input_data['p' + str(j+1)] = (donnee_bus[j][2] - donneep)/base
            input_data['q' + str(j+1)] = (donnee_bus[j][3] - donneeq)/base
            pq.append(j+1)
            v += 1
            s += 1    
        elif i[1] == 2:
            unknows.append('sigma' + str(j+1))
            unknows.append('q' + str(j+1))
            donneep = 0
            donneeq = 0
            lieu = []
            index2 = 0
            donnee_genu = []
            for i in range(len(donnee_gen)):
                donnee_genu.append(donnee_gen[i][0])
            for valeur in donnee_genu:
                if valeur == j+1:
                    donneep += donnee_gen[index2][1]
                    lieu.append(index2)
                index2 += 1
            input_data['p' + str(j+1)] = (donnee_bus[j][2] - donneep)/base
            input_data['v' + str(j+1)] = donnee_gen[compteur][5]
            slack.append(j+1)
            s += 1
            compteur += 1
        else:
            unknows.append('q' + str(j+1))
            unknows.append('p' + str(j+1))
            input_data['sigma' + str(j+1)] = 0
            input_data['v' + str(j+1)] = donnee_gen[compteur][5]
            compteur += 1
            ref.append(j+1)
        j += 1
    return slack, pq, ref, v, s, unknows, input_data

def adm(baseMVA, bus, branch):
    BUS_I, GS, BS = 0, 4, 5
    F_BUS, T_BUS, BR_R, BR_X, BR_B, TAP, SHIFT, BR_STATUS = 0, 1, 2, 3, 4,8, 9, 10 
    
    
    nb = bus.shape[0]         
    nl = branch.shape[0]       

    
    if any(bus[:, BUS_I] != list(range(1, nb+1))):
        stderr.write('buses must appear in order by bus number\n')
    stat = branch[:, BR_STATUS]              ## ones at in-service branches
    Ys = stat / (branch[:, BR_R] + 1j * branch[:, BR_X])  ## series admittance
    Bc = stat * branch[:, BR_B]              ## line charging susceptance
    tap = ones(nl)                           ## default tap ratio = 1
    i = nonzero(branch[:, TAP])              ## indices of non-zero tap ratios
    tap[i] = branch[i, TAP]                  ## assign non-zero tap ratios
    tap = tap * exp(1j * pi / 180 * branch[:, SHIFT]) ## add phase shifters
    
    Ytt = Ys + 1j * Bc / 2
    Yff = Ytt / (tap * conj(tap))
    Yft = - Ys / conj(tap)
    Ytf = - Ys / tap
    
    
    Ysh = (bus[:, GS] + 1j * bus[:, BS]) / baseMVA

    
    f = branch[:, F_BUS]                           ## list of "from" buses
    t = branch[:, T_BUS]                           ## list of "to" buses
    
    
    Cf = csr_matrix((ones(nl), (range(nl), f-ones(nl))), (nl, nb))
    
    Ct = csr_matrix((ones(nl), (range(nl), t-ones(nl))), (nl, nb))
    
    
    i = r_[range(nl), range(nl)]                   ## double set of row indices

    Yf = csr_matrix((r_[Yff, Yft], (i, r_[f-ones(nl), t-ones(nl)])), (nl, nb))
    Yt = csr_matrix((r_[Ytf, Ytt], (i, r_[f-ones(nl), t-ones(nl)])), (nl, nb))
    Ybus = Cf.T * Yf + Ct.T * Yt + \
        csr_matrix((Ysh, (range(nb), range(nb))), (nb, nb))

    return Ybus



slack, pq, ref, v, s, unknows, input_data = Ybus(donnee_bus, donnee_gen, donnee_br, nbr)
admittance = adm(base, donnee_bus, donnee_br).toarray()



start = time.time()
matrs = 2*(nbr-1) - len(slack)  #taille initial de la matrice de jordan
for e in range(1, nbr+1):
    if e not in ref:
        if 'v' + str(e) not in input_data.keys():
            globals()['x' + str(e)] = sp.symbols('x' + str(e))
        if 'sigma' + str(e) not in input_data.keys():
            globals()['z' + str(e)] = sp.symbols('z' + str(e)) 
            


diff1, diff1v = puissance(nbr, input_data, admittance, ref, pq, slack)



matricel =  [] 
for matr in range(matrs):
    matricel.append([0]*(matrs))
ange = []
for matr in range(matrs):
    ange.append([])






def calculdon(input_data, liste, x0, admittance, nbr):
    j = 2
    for i in range(nbr-1):
        input_data['sigma' + str(j)] = x0[i]
        j += 1
    i = nbr - 1
    a = 2
    while i < len(x0):
        if 'v' + str(a) in input_data:
            a += 1
        else:
            input_data['v' + str(a)] = x0[i]
            i += 1
            a += 1
    input_data['v0'] = []
    for i in range(1, nbr+1):
        input_data['v0'].append(cmath.rect(input_data['v' + str(i)], input_data['sigma' +str(i)]))
    input_data['i0'] = np.dot(admittance, input_data['v0'])
    input_data['s0'] = []
    for i in range(nbr):
        input_data['s0'].append(input_data['v0'][i] * input_data['i0'][i].conjugate())
    return input_data

def jac(diff1, nbr, ref, slack, matrice, matrs, diff1v, var, ange):
    
    c = 0
    for difft in diff1:
       
        for diffa in diff1v[c]:
            
            if str(type(diffa)) != "<class 'numpy.float64'>" and str(type(diffa)) != "<class 'numpy.int32'>" and type(diffa) != int:
                dfx = sp.diff(difft, diffa)
                if diffa in var[0]:
                    
                    matrice[c][var[1][var[0].index(diffa)]] = dfx
                    ange[c] += [var[1][var[0].index(diffa)]]
                else:
                   
                    matrice[c][var[3][var[2].index(diffa)]] = dfx
                    ange[c] += [var[3][var[2].index(diffa)]]
        c += 1
    return matrice



def detect():
    var = [[], [], [], []]
    compteur = 0
    for i in range(1, nbr+1):
        if i not in ref:
            var[0].append(globals()['z' + str(i)])
            var[1].append(compteur)
            compteur += 1
    for i in range(1, nbr+1):
        if i not in slack and i not in ref:
            var[2].append(globals()['x' + str(i)])
            var[3].append(compteur)
            compteur += 1
    return var



matrice = jac(diff1, nbr, ref, slack, matricel, matrs, diff1v, detect(), ange)



def detect2(nbr, slack):
    return [nbr, len(slack)]



def calcul(fct, var, x0):
    for i in range(len(fct)):
        if fct[i] != 0:
            fct[i] = eval(fct[i])
    return fct



def calcul1(fct, var, x0, ange, compteur):
    for i in ange[compteur]:
        fct[i] = eval(fct[i])
    return fct



def matriceobt(vecteur, ange):
    x0 = vecteur
    fx = []
    for i in diff2:
        fx.append(i) 
    fx = calcul(fx, var, x0)
    jx = [0] * matrs
    compteur = 0
    for element in range(len(matrice2)):
        jxm = []
        for i in matrice2[element]:
            jxm.append(i)
        jxm = calcul1(jxm, var, x0, ange, compteur)
        jx[element] = jxm
        compteur += 1
    return str(jx), str(fx)



def listedeliste(string):
    string = string.replace("'0'", "0")
    a = list(string)
    del a[0]
    del a[0]
    del a[-1]
    del a[-1]
    a = "".join(a)
    a = a.split('], [') 
    for i in range(len(a)):
        a[i] = a[i].split(', ')     
    for j in range(len(a)):
        for k in range(len(a)):
            if a[j][k] == "'0'":
                a[j][k] = 0
            else:
                a[j][k]  = float(a[j][k])
    return a  

    

def listestr1(string, ange):
    compteur = 0
    for j in range(len(string)):
        for i in ange[j]:
            string[j][i] = str(string[j][i])
            compteur += 1
    print(compteur)
    return string



def listestr2(string):
    a = list(string)
    del a[0]
    del a[-1]
    a = "".join(a)
    a = a.split(', ')   
    return a 


        
def liste(string):
    a = list(string)
    del a[0]
    del a[-1]
    a = "".join(a)
    a = a.split(', ')   
    for j in range(len(a)):
        a[j]  = float(a[j])
    return a



start = time.time()
x0 = np.array([0] * s + [1] * v)
sagesse = detect2(nbr, slack)
itera = 0
err = 1
errdf = 10**(-4)
var = detect()
a = str(diff1)
start = time.time()
diff2 = listestr2(a)
matrice2 = listestr1(matrice, ange)



while itera < 10 :
    if err < errdf:
        break
    for j in var[0]:
        globals()[str(var[0][var[0].index(j)])] = x0[var[1][var[0].index(j)]] 
    for j in var[2]:
        globals()[str(var[2][var[2].index(j)])] = x0[var[3][var[2].index(j)]]
    jxe, fxe = matriceobt(x0, ange)
    jxe = listedeliste(jxe)  
    fxe = liste(fxe)
    jxe = np.array(jxe)
    fxe = np.array(fxe)
    x1 = [i for i in x0]
    x0 = x0 - np.dot(np.linalg.inv(jxe), fxe)
    err = abs(max(x0 - x1))
    itera += 1
    


end = time.time()



print("calculé en ", round(end - start1, 2), "sécondes avec ", itera - 1, "itérations en newton raphson")
input_data = calculdon(input_data, sagesse, x0, admittance, nbr)
perte = np.round(sum(input_data['s0'])*base,3) 



print('='*20, 'Bus data','='*20)
print('='*50)
don = ['Bus     ', 'V       ', 'Ang     ', 's      ']
 

   
for i in don:
    if i == don[len(don)-1]:
        print(i)
    else: 
        print(i, end=' ')



for i in range(nbr):
    a = round(input_data['v' + str(i+1)], 3)
    b = round(input_data['sigma' + str(i+1)]*180/pi, 3)
    c = np.round(input_data['s0'][i]*base, 3)
    print(i+1,' '*6, a, ' '*(6-len(str(a))), b, ' '*(6-len(str(b))), c  )



print('Avec ', perte, 'comme perte de ligne du réseau')
print("calculé en ", round(end - start1, 2), "sécondes avec ", itera - 1, "itérations en newton raphson")
    

