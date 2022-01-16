# -*- coding: utf-8 -*-
import numpy as np
import pylab as plb
import csv


T=1. # Time horizon

## Loading data from vitfort.csv file
data = csv.reader(open('vitfaible.csv', newline=''), delimiter=",", quotechar='|')

## Initializing arrays
erreul_raw, liceul_raw = [], []
conterreul_raw, contliceul_raw = [], []
errmil_raw, licmil_raw = [], []
conterrmil_raw, contlicmil_raw = [], []
Npas_raw = []

## Parsing data
for row in data:
    erreul_raw.append(row[0])
    liceul_raw.append(row[1])
    conterreul_raw.append(row[2])
    contliceul_raw.append(row[3])
    errmil_raw.append(row[4])
    licmil_raw.append(row[5])
    conterrmil_raw.append(row[6])
    contlicmil_raw.append(row[7])
    Npas_raw.append(row[8])

## Data format
erreul = np.array([float(x) for x in erreul_raw[1:]])
liceul = np.array([float(x) for x in liceul_raw[1:]])
conterreul = np.array([float(x) for x in conterreul_raw[1:]])
contliceul = np.array([float(x) for x in contliceul_raw[1:]])
errmil = np.array([float(x) for x in errmil_raw[1:]])
licmil = np.array([float(x) for x in licmil_raw[1:]])
conterrmil = np.array([float(x) for x in conterrmil_raw[1:]])
contlicmil = np.array([float(x) for x in contlicmil_raw[1:]])
Npas = np.array([float(x) for x in Npas_raw[1:]])

## Pylab plots
print("Weak error Euler")
print(erreul)
print("LIC Euler")
print(liceul)
print("Weak error Euler controle variates")
print(conterreul)
print("LIC Euler controle variates")
print(contliceul)

#Plot of (h=T/N,conterreul)
plb.plot(T/Npas,conterreul, color="r", label="weak error euler")
plb.plot(T/Npas,conterreul-contliceul, color="b", label="IC 95%")
plb.plot(T/Npas,conterreul+contliceul, color="b")
plb.xlabel('Discretisation step')
plb.legend(loc="best")
plb.show()   

print("Weak error Milstein")
print(errmil)
print("LIC Milstein")
print(licmil)
print("Weak error Milstein controle variates")
print(conterrmil)
print("LIC Milstein controle variates")
print(contlicmil)

#Plot of (h=T/N,conterrmil)
plb.plot(T/Npas,conterrmil, color="r", label="weak error Milstein")
plb.plot(T/Npas,conterrmil-contlicmil, color="b", label="IC 95%")
plb.plot(T/Npas,conterrmil+contlicmil, color="b")
plb.xlabel('Discretisation step')
plb.legend(loc="best")
plb.show()  
