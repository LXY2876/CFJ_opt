import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import *
from scipy.optimize import minimize
etaE=0.6
etaFan=0.3
cd_parag=0.05
def fun(x):
    rho=1.225
    v=36.9
    S=0.9
    c=0.3
    D=c*0.21*0.5*0.5
    cl=0.0713*x[1] + 0.3327+x[2]*(1.8848*x[0]+0.3905-0.4203)
    pc=x[2]*sqrt(x[0]**2*0.01/2)*(-0.5582-12.382*x[0]**2-0.15*x[1]+x[0]**2*100/2*(1-(3/7)**2+0.2))
    cd=cd_parag+x[2]*(0.0227-1.0331*x[0]**2-0.0182-0.0022*x[1])-0.0059*cl+0.0553*cl**2
    E_c=cl/(cd/etaE+pc/etaFan)
    Pd=3000
    M_duct=x[3]
    M_motor=x[2]*0.5*1.225*v**2*pc/Pd
    m0=45
    return -E_c*(1+M_motor/m0)
def con(args):
    # 约束条件 分为eq 和ineq
    # eq表示 函数结果等于0 ； ineq 表示 表达式大于等于0  
    x1min,x1max,x2min,x2max,x3min,x3max= args
    cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x1min},\
            {'type': 'ineq', 'fun': lambda x: -x[0] + x1max},\
            {'type': 'ineq', 'fun': lambda x: x[1] - x2min},\
            {'type': 'ineq', 'fun': lambda x: -x[1] + x2max},\
            {'type': 'ineq', 'fun': lambda x: x[2] - x3min},\
            {'type': 'ineq', 'fun': lambda x: -x[2] + x3max}
                )
    return cons
 

file_path = './cal2d/统计.xlsx'

# 读取Excel文件
df = pd.read_excel(file_path)
Cmiu=df['Cmiu']
CL=df['CL']
CD=df['CD']+cd_parag
deltaCp=df['deltaCp']

# 打印数据框内容
Pc_real=df['Pc']
R=[]
Pc=[]
 
for i in range(len(Cmiu)):
    Pc.append(sqrt(Cmiu[i]*0.01/2)*(deltaCp[i]+(1-(3/7)**2+0.2)*Cmiu[i]/(2*0.01)))
     
    R.append(CL[i]/(CD[i]/etaE+Pc[i]/etaFan))
CD=np.array(CD)
Pc=np.array(Pc)
plt.figure(1)
plt.scatter(Cmiu,R)
plt.figure(2)
plt.scatter(Cmiu,CD/etaE+Pc/etaFan )
# plt.figure(3)
# plt.scatter(Cmiu,Pc_real)
print(np.polyfit(Cmiu,CD,1))
plt.figure(3)
X= np.linspace(0,0.2,100)
est=[(cd_parag+0.0227-1.0331*i)/etaE+sqrt(i*0.01/2)*(-0.5582-12.382*i+i*100/2*(1-(3/7)**2+0.2))/etaFan for i in X]
# est=[sqrt(i*0.01/2)*(-0.5582-12.382*i+i*100/2*(1-(3/7)**2)) for i in X]
plt.plot(X,est,c='r',label='est')
plt.plot(X,[sqrt(i*0.01/2)*(-0.5582-12.382*i+i*100/2*(1-(3/7)**2+0.2))/etaFan for i in X],c='b')
# plt.legend()
plt.xlabel('Cmiu')
plt.ylabel('CD+Pc')
# plt.xlim([0,0.08])
# plt.ylim([-0.02,0.02])



X_range = np.linspace(0,0.6, 100)
Y_range = np.linspace(-1, 15, 100)
X, Y = np.meshgrid(X_range, Y_range)
 
# 创建一个填充为零的 100x100 矩阵
Z = np.zeros((100, 100))

for i in range(0,100):
    for j in range(0,100):
        cl=0.0713*Y[j][1] + 0.3327+2/3*(1.8848*X[1][i]+0.3905-0.4203)
        cd=cd_parag+0.0227-1.0331*X[1][i]**2-0.0182-0.0022*Y[j][1]-0.0059*cl+0.0553*cl**2
        pc=sqrt(X[1][i]**2*0.01/2)*(-0.5582-12.382*X[1][i]**2-0.15*Y[j][1]+X[1][i]**2*100/2*(1-(3/7)**2+0.2))
        Z[j][i]= (cd/etaE+pc/etaFan)/cl
# 使用 Matplotlib 绘制网格图
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 绘制高度图
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
ax.set_xlabel('sqrt(Cmiu)')
ax.set_ylabel('alpha')
ax.set_zlabel('(cl/cd)c') 
# 添加颜色条
fig.colorbar(surf)
fig.savefig('test.tiff', dpi=1000)


cons = con((0.01,0.6,-2,6,0,0.8))
x0 = np.asarray((0.02,0,0.3))
res = minimize(fun, x0, method='SLSQP',constraints=cons)
# print(res)

plt.show()
