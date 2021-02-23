from ROOT import TFile, gROOT, TH1D, TCut
import time
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math as math
from parsing import parsing
from operator import sub
from array import array
from math import sqrt
from numpy.linalg import inv

tau = 10000.
dzeta = 0.0001
N = 188
MAT = np.eye(N)
MATC = np.eye(N)



list = [0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75, 0.7525, 0.755, 0.7575, 0.76, 0.7625, 0.765, 0.7675, 0.77, 0.7725, 0.775, 0.7775, 0.78, 0.7825, 0.785, 0.7875, 0.79, 0.7925, 0.795, 0.7975, 0.8, 0.8025, 0.805, 0.8075, 0.81, 0.8125, 0.815, 0.8175, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.9825, 0.985, 0.9875, 0.99, 0.9925, 0.995, 0.9975, 1.0, 1.0025, 1.005, 1.0075, 1.01, 1.0125, 1.015, 1.0175, 1.02, 1.0225, 1.025, 1.0275, 1.03, 1.0325, 1.035, 1.0375, 1.04, 1.0425, 1.045, 1.0475, 1.05, 1.055, 1.06, 1.065, 1.07, 1.075, 1.08, 1.09, 1.1, 1.125, 1.150, 1.175, 1.200, 1.225, 1.250, 1.275, 1.300, 1.325, 1.350, 1.375, 1.400, 1.425, 1.450, 1.475, 1.500, 1.525, 1.550, 1.575, 1.600, 1.625, 1.650, 1.675, 1.700, 1.725, 1.750, 1.775, 1.800, 1.825, 1.850, 1.875, 1.900, 1.925, 1.950, 1.975, 2.000, 2.025, 2.050, 2.075, 2.100, 2.125, 2.150, 2.175, 2.200, 2.225, 2.250, 2.275, 2.300, 2.325, 2.350, 2.375, 2.400, 2.425, 2.450, 2.475, 2.500, 2.525, 2.550, 2.575, 2.600, 2.625, 2.650, 2.675, 2.700, 2.800, 2.900, 3.000, 3.100, 3.200, 3.300, 3.400, 3.500]


def def_MATC():
    i = 0
    j = 0
    global MATC
    while i < N:
        j = 0
        while j < N:
            if i == j:
                MATC[i][j] = -1 + dzeta
                if i != 0 and i != (N-1):
                    MATC[i][j] = -2 + dzeta
            if abs(i-j)==1:
                MATC[i][j]=1
            j = j + 1;
        i = i + 1
#    MATC = inv(MATC)


def compare_thist(h1, h2):
    h1.Draw()
    h2.Draw("same")


def main(m1,m2):
    f1 = TFile("histograms/mc.root");
    h1 = TH1D("h1","h1",40,-0.04,0.04)
    h2 = TH1D("h2","h2",40,-0.04,0.04)
    ntuple = gROOT.FindObject('Tree')
    cut1 = TCut("abs(m3pimc - " + str(m1) + ") < 0.05");
    ntuple.Draw("m3pimc-m3pi >> h1",cut1)
    cut2 = TCut("abs(m3pimc - " + str(m2) + ") < 0.05");
    ntuple.Draw("m3pimc-m3pi >> h2",cut2)
    h2.Scale(h1.GetEntries()/h2.GetEntries())
    h2.SetLineColor(2)
    h1.Draw()
    h2.Draw("same")
    input()


def def_matrix():
    edata_events = experimental_events()[1]
    f1 = TFile("histograms/mc.root");
    ntuple = gROOT.FindObject('Tree')
    i = 0;
    j = 0;
    NtrueMC = array('d')
    while i < N:
        cut1 = "m3pimc > " + str(list[i]) + " && m3pimc < " + str(list[i+1]);
        Nev0 = float(ntuple.GetEntries(cut1))
        NtrueMC.append(Nev0)
        while j < N:
            if abs(i-j) < 14:
                cut2 = "m3pi > " + str(list[j]) + " && m3pi < " + str(list[j+1]);
                Nev = ntuple.GetEntries(cut2+ " && " + cut1)
                MAT[i][j] = Nev/edata_events[i] #/Nev0
            j = j+1
        i = i +1
        j = 0
    print(MAT)
    #plt.imshow(MAT)
    #plt.colorbar()
    #plt.show()
    np.savetxt('Matrix2.txt',MAT,fmt='%.1f')
    np.savetxt('NtrueMC.txt',NtrueMC,fmt='%.1f')
    input("wait")
    return 1

def events_TOY0():
    matrix = np.loadtxt('histograms/TOY.dat')
    Ndata = [0]*N
    dNdata = [0]*N
    Ndata_true = [0]*N
    i = 0
    while i < len(matrix):
        j = 0
        for el in list:
            if el < matrix[i][1]:
                j = j + 1
            else:
                break
        if j > 187:
            i = i + 1
            continue
        Ndata[j] = Ndata[j] + 1
        j = 0
        for el in list:
            if el < matrix[i][0]:
                j = j + 1
            else:
                break
        if j > 187:
            i = i + 1
            continue
        Ndata_true[j] = Ndata_true[j] + 1 
        i = i + 1
    i = 0
    for el in Ndata:
        Ndata[i] = sqrt(Ndata[i])
        dNdata[i] = Ndata[i]
        i = i + 1
    
    np.savetxt('histograms/Ndata_TOY.txt',Ndata,fmt='%.1f')
    np.savetxt('histograms/dNdata_TOY.txt',dNdata,fmt='%.1f')
    np.savetxt('histograms/Ndata_true_TOY.txt',Ndata_true,fmt='%.1f')
    return (Ndata,dNdata,Ndata_true)

def events_TOY():
    Ndata = np.loadtxt('histograms/Ndata_TOY.txt')
    dNdata = np.loadtxt('histograms/dNdata_TOY.txt')
    Ndata_true = np.loadtxt('histograms/Ndata_true_TOY.txt')
    Ndata_fold = [0]*N
    i = 0
    for el in Ndata:
        Ndata_fold[i] = Ndata[i]**2
        i = i + 1
    return (Ndata,dNdata,Ndata_true,Ndata_fold)

def experimental_events():
   # return events_TOY()

    datalow = array('d',map(sub,map(sub,map(sub,parsing("data ev/"),parsing("data e2pi/")),parsing("data e4pi/")),parsing("data bkg/")))
    ddatalow = array('d')
    i = 0
    for symb in parsing("data ev/"):
        ddatalow.append(float(sqrt(symb + parsing("data de2pi/")[i]**2 + parsing("data de4pi/")[i]**2 + parsing("data dbkg/")[i]**2)))
        i = i + 1
    Ndata = datalow + parsing("data evh/")
    dNdata = ddatalow + parsing("data devh/")
    #print(Ndata)
    i = 0
    Ndata_true = array('d')
    for el in Ndata:
        Ndata_true.append(el)
        Ndata[i] = Ndata[i]/dNdata[i]
        i = i + 1
    #print(Ndata)
    print("number of points = ",len(dNdata))
    with open('data/Ndata.txt', 'w') as file:
        i = 0
        for el in Ndata:
            file.write(str(Ndata_true[i]))
            file.write(' ')
            file.write(str(dNdata[i]))
            file.write("\n")
            i = i + 1
    return (Ndata,dNdata,Ndata_true)



def convert_to_float(MAT):
    out_array = array('d')
    for e in MAT:
     #   for e in l:
        out_array.append(np.real(e))
    out_array1 = np.asarray(out_array)
    #out_array = np.real(MAT)
    return out_array1


#    check = convert_to_float(check)
#    check = check.view(np.float64)

#    plt.imshow(check,vmax = 100)
#    print(np.real(V))
#    print(np.real(V).sort())    
#    plt.plot(np.real(V))
#    plt.colorbar()
#    plt.show()
#    input()


def cal():
    global MATC
    def_MATC()
    matrix = np.loadtxt('Matrix2.txt')
    matrix = matrix.dot(inv(MATC))
    V, W = LA.eig(matrix)
    data_events = experimental_events()[0]
    edata_events = experimental_events()[1]
    data_events_d = inv(W).dot(data_events)
    dunf_data = np.eye(N);Dunf_data = np.eye(N);
    diag_err = np.empty(N);
    energy = np.empty(N);
    i = 0
    for el in energy:
        energy[i] = list[i]/2.+list[i+1]/2.
        i = i + 1
    i = 0
    for el in data_events_d:
        data_events[i] = el*V[i]/(V[i]**2+tau)
        dunf_data[i][i] = V[i]**2/((V[i]**2+tau)**2)
        i = i + 1
    data_events = W.dot(data_events)
    data_events = inv(MATC).dot(data_events)
    ntruemc = np.loadtxt('NtrueMC.txt')
    i = 0
    for el in data_events:
        data_events[i] = el*ntruemc[i]
        i = i + 1

    dunf_data = (((inv(MATC).dot(W)).dot(dunf_data)).dot(W.transpose())).dot(inv(MATC).transpose())
    i = 0; j = 0;
    while i < N:
        j = 0
        while j < N:
            Dunf_data[i][j] = ntruemc[i]*np.real(dunf_data[i][j])*ntruemc[j]
            j = j + 1;
       # diag_err[i] = sqrt(Dunf_data[i][i])
        i = i + 1;
    #plt.plot(experimental_events()[2])
    #plt.plot(data_events)

    unf_data = convert_to_float(V) #data_events_d)
    i = 0
    for el in unf_data:
        unf_data[i] = math.fabs(el) #math.log(math.fabs(el))
        i = i + 1
    unf_data = np.sort(unf_data)
    #plt.plot(unf_data)

    #plt.imshow(Dunf_data)#,vmax = 1.2,vmin = 0.5)
    #plt.colorbar()

    plt.errorbar(energy,experimental_events()[2], yerr=edata_events, color='r',markersize = 30)
    plt.errorbar(energy,data_events, yerr=diag_err, color='b',markersize = 30)
    #plt.errorbar(energy,experimental_events()[3], color='g',markersize = 30)
    plt.show()
    return 1



#print(experimental_events()[0])
#events_TOY()
#cal()
experimental_events()
#main(0.76,1.72)
#def_matrix()
#def_MATC()
#print(inv(MATC))

