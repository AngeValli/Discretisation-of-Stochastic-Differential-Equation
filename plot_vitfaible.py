import numpy as np
import pylab as plb
import csv


T=1. #horizon

data = csv.reader(open('vitfaible.csv', newline=''), delimiter=",", quotechar='|')
erreul_raw, liceul_raw = [], []
conterreul_raw, contliceul_raw = [], []
errmil_raw, licmil_raw = [], []
conterrmil_raw, contlicmil_raw = [], []
Npas_raw = []

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

erreul = np.array([float(x) for x in erreul_raw[1:]])
liceul = np.array([float(x) for x in liceul_raw[1:]])
conterreul = np.array([float(x) for x in conterreul_raw[1:]])
contliceul = np.array([float(x) for x in contliceul_raw[1:]])
errmil = np.array([float(x) for x in errmil_raw[1:]])
licmil = np.array([float(x) for x in licmil_raw[1:]])
conterrmil = np.array([float(x) for x in conterrmil_raw[1:]])
contlicmil = np.array([float(x) for x in contlicmil_raw[1:]])
Npas = np.array([float(x) for x in Npas_raw[1:]])

print("Erreur faible Euler")
print(erreul)
print("LIC Euler")
print(liceul)
print("Erreur faible Euler VA controle")
print(conterreul)
print("LIC Euler VA controle")
print(contliceul)

#representation graphique de (h=T/N,conterreul)
plb.plot(T/Npas,conterreul, color="r", label="erreur faible euler")
plb.plot(T/Npas,conterreul-contliceul, color="b", label="IC 95%")
plb.plot(T/Npas,conterreul+contliceul, color="b")
plb.xlabel('Pas de discretisation')
plb.legend(loc="best")
plb.show()   

print("Erreur faible Milstein")
print(errmil)
print("LIC Milstein")
print(licmil)
print("Erreur faible Milstein VA controle")
print(conterrmil)
print("LIC Milstein VA controle")
print(contlicmil)

#representation graphique de (h=T/N,conterreul)
plb.plot(T/Npas,conterrmil, color="r", label="erreur faible Milstein")
plb.plot(T/Npas,conterrmil-contlicmil, color="b", label="IC 95%")
plb.plot(T/Npas,conterrmil+contlicmil, color="b")
plb.xlabel('Pas de discretisation')
plb.legend(loc="best")
plb.show()  
