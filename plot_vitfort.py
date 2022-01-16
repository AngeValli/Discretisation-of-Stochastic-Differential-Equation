# -*- coding: utf-8 -*-
import numpy as np
import pylab as plb
import csv


T=1. # Time horizon

## Loading data from vitfort.csv file
data = csv.reader(open('vitfort.csv', newline=''), delimiter=",", quotechar='|') 

## Initializing arrays
erreul_raw, liceul_raw = [], []
errmil_raw, licmil_raw = [], []
Npas_raw = []

## Parsing data
for row in data:
    erreul_raw.append(row[0])
    liceul_raw.append(row[1])
    errmil_raw.append(row[2])
    licmil_raw.append(row[3])
    Npas_raw.append(row[4])

## Data format
erreul = np.array([float(x) for x in erreul_raw[1:]])
liceul = np.array([float(x) for x in liceul_raw[1:]])
errmil = np.array([float(x) for x in errmil_raw[1:]])
licmil = np.array([float(x) for x in licmil_raw[1:]])
Npas = np.array([float(x) for x in Npas_raw[1:]])

## Pylab plots
plb.plot(T/Npas,erreul, color="r", label="Quadratic error Euler")
plb.plot(T/Npas,erreul-liceul, color="b", label="IC 95%")
plb.plot(T/Npas,erreul+liceul, color="b")
plb.xlabel('Discretisation step')
plb.legend(loc="best")
plb.show()

#Display of error vectors and half-width of IC at 95%
print("Error euler")
print(erreul)
print("LIC Euler")
print(liceul)
plb.plot(T/Npas**2,errmil, color="r", label="Quadratic error Milstein")
plb.plot(T/Npas**2,errmil-licmil, color="b", label="IC 95%")
plb.plot(T/Npas**2,errmil+licmil, color="b")
plb.xlabel('Square of discretisation step')
plb.legend(loc="best")
plb.show()

#Display of error vectors and half-width of IC at 95%
print("Error Milstein")
print(errmil)
print("LIC Milstein")
print(licmil)
