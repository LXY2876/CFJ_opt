import numpy as np
from math import *
 
from scipy.optimize import minimize
import sys
etaE=0.6
etaFan=0.48
W_S=105
# cd_parag=0.00646
cd_parag=0.00646
rho=1.225
# v=36.8
# S=0.9
# c=0.3
v=100
S=0.9
c=0.3
g=9.8
D=c*0.21*0.5
t_duct=0.001
rho_duct=1200

Pd=3000
CL=0.620996703
E=CL/(cd_parag+0.037030164)
m0=0.5*1.225*S*CL*v**2/g
m_payload=0.2*m0
 
hinj=c*0.01
hsuc=c*7/300
SF=0.3
N=6
global global_best
def fun(x):
    if x[0]<0:
        return sys.maxsize
    cl= 0.620996703+x[1]*x[2]*D*c/S*(2.2795*sqrt(x[0])+0.6973-0.759519745)
    k=0.0259*x[2]**2-0.1178*x[2]+0.2967
    pc=x[1]*x[2]*D*c/S*sqrt(x[0]*0.01/2)*(-0.662-17.327*x[0]+x[0]*100/2*(1-(3/7)**2+k))
    # cd=cd_parag+0.037030164+x[1]*x[2]*D/S*(0.0297-0.91115*x[0]-0.0254061)-0.0059*cl+0.0553*cl**2
    cd=cd_parag+0.037030164+x[1]*x[2]*D*c/S*(0.0297-0.91115*x[0]-0.0254061)
    E_c=cl/(cd+pc*etaE/etaFan)
    M0=0.5*rho*S*v**2*cl/g
    M_duct=x[1]*((hinj+hsuc)/2+x[2]*D+np.pi*D)*c*rho_duct*t_duct
    M_motor=0.5*1.225*v**3*S*pc/Pd /etaFan
    return (14.27901503*(1-m_payload/m0-SF))/(E_c*(1-(m_payload+M_duct+M_motor)/M0-SF)),M_motor,M_duct,cl

    # return (14.27901503*(1-m_payload/m0-SF))/(E_c*(1-(m_payload)/m0-SF)) 
# def fun1(x):
#     # N=1
#     cl= 0.620996703+N*x[1]*D/S*(2.2795*sqrt(x[0])+0.6973-0.759519745)
#     k=0.0259*x[1]**2-0.1178*x[1]+0.2967
#     pc=N*x[1]*D/S*sqrt(x[0]*0.01/2)*(-0.662-17.327*x[0]+x[0]*100/2*(1-(3/7)**2+k))
#     # cd=cd_parag+0.037030164+x[1]*x[2]*D/S*(0.0297-0.91115*x[0]-0.0254061)-0.0059*cl+0.0553*cl**2
#     cd=cd_parag+0.037030164+N*x[1]*D/S*(0.0297-0.91115*x[0]-0.0254061)
#     E_c=cl/(cd/etaE+pc/etaFan)
#     M0=0.5*rho*S*v**2*cl/g
#     M_duct=N*((hinj+hsuc)/2+x[1]*D+np.pi*D)*c*rho_duct*t_duct
#     M_motor=0.5*1.225*v**3*S*pc/Pd
#     return (14.27901503*etaE*(1-m_payload/m0-SF))/E_c*(1-(m_payload+M_duct+M_motor)/M0-SF)
# def con1(args):
#     # 约束条件 分为eq 和ineq
#     # eq表示 函数结果等于0 ； ineq 表示 表达式大于等于0  
#     x1min,x2min,x1max,x2max= args
#     cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x1min},\
#             {'type': 'ineq', 'fun': lambda x: -x[0] + x1max},\
#             {'type': 'ineq', 'fun': lambda x: x[1] - x2min},\
#             {'type': 'ineq', 'fun': lambda x: -x[1] + x2max},\
#             {'type': 'ineq', 'fun': lambda x: 0.8 - N*x[1]*D*c/S}
#                 )
#     return cons
def con(args):
    # 约束条件 分为eq 和ineq
    # eq表示 函数结果等于0 ； ineq 表示 表达式大于等于0  
    x1min,x2min,x3min,x1max,x2max,x3max= args
    cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x1min},\
            {'type': 'ineq', 'fun': lambda x: -x[0] + x1max},\
            {'type': 'ineq', 'fun': lambda x: x[1] - x2min},\
            {'type': 'ineq', 'fun': lambda x: -x[1] + x2max},\
            {'type': 'ineq', 'fun': lambda x: x[2] - x3min},\
            {'type': 'ineq', 'fun': lambda x: -x[2] + x3max},\
            {'type': 'ineq', 'fun': lambda x: 0.8 - x[2]*x[1]*D*c/S}
                )
    return cons

global_best = {'val': inf, 'x': None}  # 存储最优目标函数值和对应的变量值

def integerPro(args,t=1.0E-7):
	# 求解松弛问题
    local_args=args.copy()
    # if(args[1]==6):
    print(args)
    cons = tuple(local_args)
    x0 = np.array(local_args[0:3]) + np.array(local_args[3:6])
    x0 = x0 / 2   
     
    res = minimize(fun, x0, method='SLSQP',constraints=con(cons))
	# res = linprog(c, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq)
    bestVal = res.fun
    bestX = res.x
    # if not(type(res.x) is float or res.success != 0): 
    #     bestVal = res.fun
    # else:
    #     bestVal=sys.maxsize
    # 停止条件 & bound
    if res.success==False or bestVal >= global_best['val']:
        return (global_best['val'], global_best['x'])
 

    # 当前解成为新的全局最优解
    if bestVal < global_best['val'] and (((bestX[1]-floor(bestX[1]))<t or (ceil(bestX[1])-bestX[1])<t)):
        global_best['val'] = bestVal
        global_best['x'] = bestX

    if ((bestX[1]-floor(bestX[1]))<t or (ceil(bestX[1])-bestX[1])<t) :
        return (bestVal,bestX)
    else:
    	# 进行branch，这里简单选择第一个非整数变量
        
        bound1=floor(bestX[1])
        bound2=ceil(bestX[1])
        cons1 = args.copy()
        cons1[1]=bound2
        cons2 = args.copy()
        cons2[4]=bound1

        r1 = integerPro(cons1,t)
        r2 = integerPro(cons2,t)
        # tree search，这里使用width first
        if r1[0] <r2[0]:
            return r1
        else:
            return r2
# print(integerPro(np.array([0.01,0,3.31,0.2,6,6.62])))
if __name__== "__main__" :
    # print(integerPro(np.array([0.01,0,3.31,0.2,26,6.62])))
    x0=np.array([0.05,4,4])
    # res=minimize(fun, x0, method='SLSQP',constraints=con((0.01,1,3.31,0.2,26,6.62)))
    # x0=np.array([0.05,6,4])
    # res=integerPro(np.array([0.01, 1, 3.31, 0.2, 26, 6.62]))
    # res=minimize(fun, x0, method='SLSQP',constraints=con((0.01,6,3.31,0.2,6,6.62)))
    # print(res)
    # print(fun([0.06869002,25,3.31260421]))
    
    data = [
        [0.09189827, 1, 3.66374923],
        [0.09020302, 13, 3.74911453],
        [0.08833897, 20, 3.77517097],
        [0.08668294, 20, 3.80069177],
        [0.08785143, 21, 3.62811791],
        [0.08592648, 22, 3.46320346],
        [0.08136069, 23, 3.3126294],
        [0.07281302, 23, 3.3126294],
        [0.07081615, 23, 3.3126294],
        [0.06957445, 23, 3.3126294],
        [0.06895199, 23, 3.3126294],
        [0.06744462, 23, 3.3126294]
    ]
    results = [fun(x) for x in data]

# # 打印结果
    for result in results:
        print(result)