# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 19:46:46 2023

@author: ULRICH_LUMENDO
"""

from fonction import ordre, Spuissance, dSpuissance_dV, methode, admitt, execfp

while True:
    casetype = input('Quel fichier de données voulez-vous utilser : ')
    if casetype == '':
        break
    try:
        execfp(casetype)
    except Exception:
        print(f"Le fichier {casetype}.py n'est pas disponible")
    


    

    

    
    
    
        