import sympy as sp
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
