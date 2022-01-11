# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 16:24:10 2014

@author: jourdain
"""

import numpy as np
import pylab as plb
import csv


T=1. #horizon

data = csv.reader(open('vitfort.csv', newline=''), delimiter=",", quotechar='|')
erreul_raw, liceul_raw = [], []
errmil_raw, licmil_raw = [], []
Npas_raw = []

for row in data:
    erreul_raw.append(row[0])
    liceul_raw.append(row[1])
    errmil_raw.append(row[2])
    licmil_raw.append(row[3])
    Npas_raw.append(row[4])
    

erreul = np.array([float(x) for x in erreul_raw[1:]])
liceul = np.array([float(x) for x in liceul_raw[1:]])
errmil = np.array([float(x) for x in errmil_raw[1:]])
licmil = np.array([float(x) for x in licmil_raw[1:]])
Npas = np.array([float(x) for x in Npas_raw[1:]])


plb.plot(T/Npas,erreul, color="r", label="Erreur quadratique Euler")
plb.plot(T/Npas,erreul-liceul, color="b", label="IC 95%")
plb.plot(T/Npas,erreul+liceul, color="b")
plb.xlabel('Pas de discretisation')
plb.legend(loc="best")
plb.show()
#Affichage des vecteurs d'erreurs et des demi-largeurs des IC à 95%
print("Erreur euler")
print(erreul)
print("LIC Euler")
print(liceul)
plb.plot(T/Npas**2,errmil, color="r", label="Erreur quadratique Milstein")
plb.plot(T/Npas**2,errmil-licmil, color="b", label="IC 95%")
plb.plot(T/Npas**2,errmil+licmil, color="b")
plb.xlabel('Carre du pas de discretisation')
plb.legend(loc="best")
plb.show()
#Affichage des vecteurs d'erreurs et des demi-largeurs des IC à 95%
print("Erreur Milstein")
print(errmil)
print("LIC Milstein")
print(licmil)
