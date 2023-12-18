from math import *
from NSGA2 import *
import numpy as np
import pandas as pd
import pickle
#X=[N,W_D,Cmiu_takeoff,Cmiu_curise]


class PLANE():
     
    def __init__(self,x):

        self.x=x

        self.rho=1.225
        self.c=0.3
        self.A=10
        self.S=self.c**2*self.A
        self.structure=0.5
        self.W_S=600
        self.G0=self.S*self.W_S#飞机总重(N)
        self.G1=0
        self.k=1/(pi*self.A*0.8) #奥斯瓦尔德效率因子取0.8
        self.parasite_drag=0.02095
        self.CL_tof0=1
        self.CL_tof=1.2
        self.V_tof=0
        self.CL_curise=0.55
        self.CD_curise=0.02
        self.V=sqrt(self.W_S/(0.5*self.rho*self.CL_curise))
        self.Cd0_tof=0.1
        self.CD0_curise=0.1
        self.beta=1
        
        self.Cmiu_takeoff=x[2]
        self.Cmiu_curise=x[3]
        self.Pd=3000 #w/kg
        self.eta=0.1
        self.eta_engine=0.7
        self.T_max=0.5*self.G0
        self.Gb=0.4*self.G0
        self.Eb=250 #Wh/kg
        self.Duct=CFJ(x,self.c)
        self.S_CFJ=x[0]*self.c*self.Duct.W
        self.miu=0.2#地面滑跑摩擦阻力系数
        self.TofDistance=0
        self.Range=0
        self.tof_ref=0
        self.curise_ref=0
    def get_raw(self):
        self.V_tof=1.2*sqrt(self.G0/(0.5*self.rho*self.CL_tof))
        F0=self.T_max-self.miu*self.G0
        CD=self.Cd0_tof+self.k*self.CL_tof**2
        F1=self.T_max-0.5*self.V_tof**2*self.S*(CD)

        self.tof_ref=self.G0/(2*9.8)*self.V_tof**2/(F0-F1)*log(F0/F1)
        CD_curise=self.CD0_curise+self.k*self.CL_curise**2
        self.curise_ref=self.Gb*self.Eb/(self.G0)*(self.CL_curise/(CD_curise/self.eta_engine )) 
        # return [tof ,R]
    def update_aero(self):
        self.get_raw()
        self.CL_tof=self.CL_tof0+6*self.S_CFJ/self.S*self.Cmiu_takeoff
        self.CL_curise=self.CL_curise+6*self.S_CFJ/self.S*self.Cmiu_curise
        self.CD_curise=self.CD0_curise-self.beta*self.S_CFJ/self.S*self.Cmiu_curise+self.k*self.CL_curise**2
         
    def update_G(self):
        try:
            # a=cardano_solution(0.5*self.rho*self.S_CFJ*self.Duct.pc_tof/(self.Pd*self.eta)*9.8,-0.5*self.rho*self.S*self.CL_tof,0,self.G0+self.Duct.duct)
            roots = np.roots([0.5*self.rho*self.S_CFJ*self.Duct.pc_tof/(self.Pd*self.eta)*9.8,-0.5*self.rho*self.S*self.CL_tof, 0, self.G0+self.Duct.duct])

            # 选出实数根
            real_roots = [np.real(root) for root in roots if np.isreal(root)]
        except:
            print("求根失败",self.x )
            exit()  
        # print(real_roots)
        
        self.V_tof=min([i for i in real_roots if i>0])
        motor_curise=0.5*self.rho*self.V**3*self.S*self.Duct.pc_curise*9.8/(self.Pd*self.eta)
        motor_tof=0.5*self.rho*self.V_tof**3*self.S_CFJ*self.Duct.pc_tof*9.8/(self.Pd*self.eta)
        # print("duct",self.Duct.duct)
        # print("motor_curise",motor_curise)
        # print("motor_tof",0.5*self.rho*self.S_CFJ*self.Duct.pc_tof/(self.Pd*self.eta)*self.V_tof**3*9.8)
        self.G1=max(motor_tof,motor_curise)+self.G0+self.Duct.duct
    


    def get_performance(self):
         
        F0=self.T_max-self.miu*self.G1
        CD=self.Cd0_tof-self.beta*self.Cmiu_takeoff+self.k*self.CL_tof**2
        F1=self.T_max-0.5*self.V_tof**2*self.S*(CD)
        # print(F0,F1)
        try:
            self.TofDistance=self.G1/(2*9.8)*self.V_tof**2/(F0-F1)*log(F0/F1)
        except:
            print("起飞距离计算失败",F0,F1,self.G1,self.x)
            self.TofDistance=1000
            # exit()
        # CD_curise=self.CD0_curise-self.beta*self.Cmiu_curise
        self.Range=self.Gb*self.Eb/(self.G1)*(self.CL_curise/(self.CD_curise/self.eta_engine+self.S_CFJ/self.S*self.Duct.pc_curise/self.eta)) 
        return [self.TofDistance,self.Range]
    
    def update_date(self,x):
        
        pass

class CFJ():
    def  __init__(self,x,c):
         
        self.rho=1.225
        self.D=0.1*c
        self.c=c #弦长
        self.N=x[0]
        self.W=x[1]*self.D
        self.Cmiu_tof=x[2]
        self.Cmiu_curise=x[3]
        self.deltaCP0_Tof=-1
        self.deltaCP0_curise=-1
        self.deltaCP_K=-20
        self.rho_duct=1050

        self.t=0.001#mm
        self.inj_c=0.01
        self.suc_c=7/300
        # with open('../duct_design/gpr_model.pkl', 'rb') as file:
        #     self.model_gpr = pickle.load(file)
        # self.K1,_=self.model_gpr.predict([[x[1],x[2]]], return_std=True)#管道总压损失系数
        # self.K2,_=self.model_gpr.predict([[x[1],x[3]]], return_std=True)#管道总压损失系数
        # self.K1=self.K1[0]
        # self.K2=self.K2[0]
        self.K1=0.0*(self.W/self.D*0.1+0.018*(self.W/self.D)**2)#管道总压损失系数
        self.K2=0.0*(self.W/self.D*0.1+0.018*(self.W/self.D)**2)#管道总压损失系数
        self.duct=x[0]*self.duct_mass()*9.8#管道重量(N)
        self.pc_tof=self.Pc_Takeoff()
        self.pc_curise=self.Pc_Curise()

    def duct_mass(self):
        return self.rho_duct*(2*self.W+self.D*pi)*self.c*self.t
    def Pc_Takeoff(self):
        # print(self.K1)
        return  sqrt(self.Cmiu_tof*self.inj_c/2)*\
            (self.deltaCP0_Tof+self.deltaCP_K*self.Cmiu_tof+self.Cmiu_tof/(2*self.inj_c))*(1-(self.inj_c/self.suc_c)**2+self.K1)
    def Pc_Curise(self):
        return sqrt(self.Cmiu_curise*self.inj_c/2)*\
            (self.deltaCP0_curise+self.deltaCP_K*self.Cmiu_curise+self.Cmiu_curise/(2*self.inj_c))*(1-(self.inj_c/self.suc_c)**2+self.K2)

 

def function1(x):
    plane=PLANE(x)
    plane.update_aero()
    plane.update_G()
    
    plane.get_performance()
    return  1/(plane.TofDistance /plane.tof_ref)

def function2(x):
    plane=PLANE(x)
    plane.update_aero()
    plane.update_G()
    plane.get_performance()
    return  plane.Range / plane.curise_ref
def CV_value(x):
    # plane=PLANE(x)
    # if x[0]*plane.Duct.W<=0.8*plane.A*plane.c+0.001*plane.A*plane.c:
    
    #     return 0
    # else:
    #     return (x[0]*plane.Duct.W-0.8*plane.A*plane.c)
    return 0

if __name__=='__main__':
    nsga=NSGA2(function1,function2,CV_value,[1,0,0,0])
    nsga.bounds=[(2,14),(2,7),(0.08,0.3),(0.03,0.1)]
    result,F1,F2=nsga.main()
    v1=max(F1)
    v2=max(F2)
    baseline=PLANE([0,0,0,0]).get_raw()
    # v1=baseline[0]
    # v2=-1*baseline[1]
    v1=-1
    v2=-1
    index=0
    max_obj=0
    
    # F1=np.array(F1)/baseline[0]
    # F2=np.array(F2)/baseline[1]
    for i in range(0,nsga.pop_size):
        if (F1[i]-v1)*(F2[i]-v2)>max_obj:
            max_obj=(F1[i]-v1)*(F2[i]-v2)
            index=i
    X=np.linspace(F1[index]-0.2,F1[index]+0.2,30)
    Y=max_obj/(X-v1)+v2

    planes=[PLANE(x) for x in result]
    other_info=np.zeros([nsga.pop_size,3])
    i=0
    for plane in planes:
        plane.update_aero()
        plane.update_G()
        plane.get_performance()
        other_info[i][0]=plane.G1
        other_info[i][1]=plane.Duct.duct
        other_info[i][2]=plane.G1-plane.Duct.duct-plane.G0
        i=i+1
    print("nash",result[index],F1[index],F2[index])
    plt.scatter(X,Y, marker='^',color='r') 
    
    plt.xlabel('Tof', fontsize=15)
    plt.ylabel('Curise', fontsize=15)
    plt.scatter(F1, F2, marker='+',color='b')
    # plt.xlim(0, 1)
    # plt.ylim(-2, -1)
    dic1={"N":result[:,0],"W_D":result[:,1],"Cmiu_Tof":result[:,2],"Cmiu_Curise":result[:,3],"TofDistance":F1,"CuriseRange":F2,"总重":other_info[:,0],"管道重量":other_info[:,1],"电机重量":other_info[:,2]}
    df=pd.DataFrame(dic1)
    df.to_excel("结果3.xlsx",sheet_name="sheet1",index=False)
    plt.show()

    # plane=PLANE([2.00037998 ,0.06   ,    0.19845694, 0.03937537])
    # print(plane.get_raw())
    # plane.update_aero()
    # plane.update_G()
    # print(plane.G1)
    # print(plane.get_performance()) 
