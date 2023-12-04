

 
#X=[N,W.Cmiu_takeoff,Cmiu_curise]
class PLANE():
    def __init__(self,x):

        self.Cmiu_takeoff=x[2]
        self.Cmiu_curise=x[3]
        self.G=100
        self.A=10
        
        self.c=0.3
        
        self.motor=0
        
        self.pd=3000


class CFJ():

    def  __init__(self,x):
        self.N=x[0]
        self.W=x[1]
        self.rho_duct=1050
        self.D=0.03
        self.c=0.3
        self.duct=self.duct_mass()
    def duct_mass(self):
        self.duct=self.rho_duct*(2*self.W+self.D*pi)*self.c
 

if __name__=='__main__':
    print(pi)
    plane=CFJ([1,2,3,4])
    print(plane.duct)
     
