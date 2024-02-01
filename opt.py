import sys
import numpy as np
sys.path.append('../')
from NSGA2 import *
from math import *
etaE=0.6
etaFan=0.45
cd_parag=0.023
rho=1.225
v=36.9
S=0.9
c=0.3
D=c*0.21*0.5
t_duct=0.001
rho_duct=1200

Pd=3000
m0=45
s=0.3
hinj=c*0.01
hsuc=c*7/300
def fun1(X):#sqrt(Cmiu),alpha,N,W/D
    pop, x_num = X.shape
    y1=[0 for i in range(pop)]
    y2=[0 for i in range(pop)]
 
    CV3=[0 for i in range(pop)]
    for i in range(pop):
        x=X[i]
        cl=0.0713*x[1] + 0.3327+x[2]*x[3]*D*(1.8848*x[0]+0.3905-0.4203)
        pc=x[2]*x[3]*D*sqrt(x[0]**2*0.01/2)*(-0.5582-12.382*x[0]**2-0.15*x[1]+x[0]**2*100/2*(1-(3/7)**2+0.2))
        cd=cd_parag+x[2]*x[3]*D*(0.0227-1.0331*x[0]**2-0.0182-0.0022*x[1])-0.0059*cl+0.0553*cl**2
        E_c=cl/(cd/etaE+pc/etaFan)
        
        M_duct=x[2]*((hinj+hsuc)/2+x[3]*D+np.pi*D)*c*rho_duct*t_duct
        M_motor=0.5*1.225*v**2*pc/Pd
        y1[i]=E_c/(1+(M_motor+M_duct)/m0/(1-s))
        y2[i]=0
        if x[2]*x[3]*D/S<0.8:
            CV3[i]=0
        else:
            CV3[i]=x[2]*x[3]*D/S-0.8

    return y1, y2,CV3
def fun2(x):
    return 0
def CV_value():
    return 0
if __name__=='__main__':
      
    nsga=NSGA2(fun1,fun2,fun1,CV_value,[0,0,1,0])
    nsga.max_gen=100
    nsga.bounds=[(0.01,0.6),(-2,7),(0,100),(2,7)]
    result,F1,F2=nsga.main()
    print(result,F1)