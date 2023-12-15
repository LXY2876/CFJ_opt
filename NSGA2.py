# Program Name: NSGA-II.py
# author；lxy
# Time:2023/03/11
 
#Importing required modules
from math import *
import random
import matplotlib.pyplot as plt
import numpy as np
import sys
from numpy import pi
#First function to optimize

class NSGA2():
    def __init__(self,func1,func2,constraint,flag):
    #Initialization

        self.min_x=[-55,-55]
        self.max_x=[55,55]
        self.flag=flag
        self.pop_size=40
        self.max_gen=400
        self.obj1=func1
        self.obj2=func2
        self.constraint=constraint
        self.bounds=[(-55,55),(-55,55)]
        self.tour_size=2
    # def obj1(self,x):
    #     value = -x[0]**2
    #     return value
    
    # #Second function to optimize
    # def obj2(self,x):
    #     value = -(x[0]-2)**2
    #     return value
    
    #Function to find index of list,且是找到的第一个索引
    def first_population(self):
        p=13
        r=[]
        solution=[]
        s=len(self.bounds)
        for k in range(1,s+1):
            r.append(2*cos(2*pi*k/p))
        r=np.array(r)
        for k in range(1,self.pop_size+1):
            P=(k*r)%1
            for i in range(s):
                P[i]=P[i]*(self.bounds[i][1]-self.bounds[i][0])+self.bounds[i][0]
            solution.append(P.tolist())
        return solution
    def index_of(self,a,list):
        for i in range(0,len(list)):
            if list[i] == a:
                return i
        return -1
    
    #Function to sort by values 找出front中对应值的索引序列
    def sort_by_values(self,list1, values):
        sorted_list = []
        while(len(sorted_list)!=len(list1)):
            if self.index_of(min(values),values) in list1:
                sorted_list.append(self.index_of(min(values),values))
            values[self.index_of(min(values),values)] = inf
        return sorted_list
    
    #Function to carry out NSGA-II's fast non dominated sort
    def fast_non_dominated_sort(self,values1, values2,CV_value):
        S=[[] for i in range(0,len(values1))]   #len(values1)个空列表
        front = [[]]
        n=[0 for i in range(0,len(values1))]
        rank = [0 for i in range(0, len(values1))]
        #将front0全部整理出来了，并未对front1-n等进行整理
        for p in range(0,len(values1)):
            S[p]=[]
            n[p]=0
            for q in range(0, len(values1)):

                if CV_value[p]==0 and CV_value[q]==0:
                     
                    if (values1[p] > values1[q] and values2[p] > values2[q]) or (values1[p] >= values1[q] and values2[p] > values2[q]) or (values1[p] > values1[q] and values2[p] >= values2[q]):
                        if q not in S[p]:
                            S[p].append(q)
                    elif (values1[q] > values1[p] and values2[q] > values2[p]) or (values1[q] >= values1[p] and values2[q] > values2[p]) or (values1[q] > values1[p] and values2[q] >= values2[p]):
                        n[p] = n[p] + 1
                elif CV_value[p]==0 and CV_value[q]>0:
                    if q not in S[p]:
                        S[p].append(q)    
                elif CV_value[p]>0 and CV_value[q]==0:
                    n[p] = n[p] + 1            

                else:
                    if CV_value[p]<CV_value[q]:
                        if q not in S[p]:
                            S[p].append(q)
                    elif CV_value[p]>CV_value[q]:
                        n[p] = n[p] + 1     
            if n[p]==0:
                rank[p] = 0
                if p not in front[0]:
                    front[0].append(p)
        i = 0
        #该循环能将所有的个体全部进行分类，显然最后一层的个体中，没有可以支配的个体了
        while(front[i] != []):
            Q=[]
            for p in front[i]:
                for q in S[p]:
                    n[q] =n[q] - 1
                    if( n[q]==0):
                        rank[q]=i+1
                        if q not in Q:
                            Q.append(q)
            i = i+1
            front.append(Q)
    
        del front[len(front)-1]    #删除了最后一层无支配个体的front层,最后一层是空集
        return front
    
    #Function to calculate crowding distance  同层之间的一个计算，给出一层中的拥挤度
    def crowding_distance(self,values1, values2, front):
        # distance = [0 for i in range(len(front))]
        lenth= len(front)
        for i in range(lenth):
            distance = [0 for i in range(lenth)]
            sorted1 = self.sort_by_values(front, values1[:])  #找到front中的个体索引序列
            sorted2 = self.sort_by_values(front, values2[:])  #找到front中的个体索引序列
            distance[0] = 4444
            distance[lenth-1] = 4444
            for k in range(1,lenth-1):
                distance[k] = distance[k]+ (values1[sorted1[k+1]] - values1[sorted1[k-1]])/(max(values1)-min(values1))
                # print("/n")
                # print("k:",k)
                # print("distance[{}]".format(k),distance[k])
            for k in range(1,lenth-1):
                distance[k] = distance[k]+ (values2[sorted2[k+1]] - values2[sorted2[k-1]])/(max(values2)-min(values2))
        return distance
    
    # #Function to carry out the crossover
    def crossover(self,a,b):
        r=random.random()
        if r>0.5:
            return self.mutation((a+b)/2)
        else:
            return self.mutation((a-b)/2)
    # #Function to carry out the mutation operator
    def mutation(self,solution):
        mutation_prob = random.random()
        if mutation_prob <1:
            solution[0] = self.min_x[0]+(self.max_x[0]-self.min_x[0])*random.random()
            solution[1] = self.min_x[1]+(self.max_x[1]-self.min_x[1])*random.random()
        return solution
    # def tournament_selection(self,non_dominated_sorted_solution,crowding_distance_values):
    #     pool_size=self.pop_size/2
    #     candidate=[0 for i in range(self.tour_size)]
    #     for i in range(pool_size):
    #         for j in range(self.tour_size):
    #             candidate[j]=random.randint(0, self.pop_size - 1)
    #             if j>0:
    #                 while candidate[j] not in candidate[0:j-1]:
    #                     candidate[j]=random.randint(0, self.pop_size - 1)
    #         candidate_rank=[for ]
          

    def cross_mutation(self,chromo_parent,  pc, pm, yita1, yita2):
        # Simulated binary crossover and polynomial mutation
        pop, x_num = chromo_parent.shape
        suoyin = 0
        off = np.zeros_like(chromo_parent)

        for i in range(int(pop/2)):
            # Simulated binary crossover
            # Initialize offspring population
            # Randomly select two parent individuals
            parent_1 = random.randint(0, pop - 1)
            parent_2 = random.randint(0, pop - 1)
            
            # Ensure two parent individuals are not the same
            while np.array_equal(chromo_parent[parent_1, :], chromo_parent[parent_2, :]):
                parent_2 = random.randint(0, pop - 1)

            chromo_parent_1 = chromo_parent[parent_1, :]
            chromo_parent_2 = chromo_parent[parent_2, :]
            off_1 = chromo_parent_1.copy()
            off_2 = chromo_parent_2.copy()

            if random.random() < pc:
                # Perform simulated binary crossover
                u1 = np.random.rand(x_num)
                gama = np.zeros(x_num)
                for j in range(x_num):
                    if u1[j] < 0.5:
                        gama[j] = (2 * u1[j]) ** (1 / (yita1 + 1))
                    else:
                        gama[j] = (1 / (2 * (1 - u1[j]))) ** (1 / (yita1 + 1))

                    off_1[j] = 0.5 * ((1 + gama[j]) * chromo_parent_1[j] + (1 - gama[j]) * chromo_parent_2[j])
                    off_2[j] = 0.5 * ((1 - gama[j]) * chromo_parent_1[j] + (1 + gama[j]) * chromo_parent_2[j])
                    
                    # Ensure offspring are within defined domain
                    off_1[j] = np.clip(off_1[j], self.bounds[j][0], self.bounds[j][1])
                    off_2[j] = np.clip(off_2[j], self.bounds[j][0], self.bounds[j][1])
                # print("gama",gama)    
        

            # Polynomial mutation
            etm=20
            if random.random() < pm:
                u2 = np.random.rand(x_num)
                delta = np.zeros(x_num)
                for j in range(x_num):
                    if random.random() < 1/x_num:
                    # if u2[j] < 0.5:
                    #     delta[j] = (2 * u2[j]) ** (1 / (yita2 + 1)) - 1
                    # else:
                    #     delta[j] = 1 - (2 * (1 - u2[j])) ** (1 / (yita2 + 1))

                        u=np.random.rand()
                        delta1=(off_1[j]-self.bounds[j][0])/(self.bounds[j][1]-self.bounds[j][0])
                        delta2=(self.bounds[j][1]- off_1[j])/(self.bounds[j][1]-self.bounds[j][0])
                        if u2[j] <0.5:
                            delta[j]=pow(2*u2[j]+(1-2*u2[j])*(1-delta1)**(etm+1),1/(etm+1))-1
                        else:
                            delta[j]=1-pow(2*(1-u2[j])+2*(u2[j]-0.5)*(1-delta2)**(etm+1),1/(etm+1))

                        off_1[j] += delta[j]*(self.bounds[j][1]-self.bounds[j][0])
                        # print("delta",delta[j])
                        # off_1[j] += delta[j]
                        off_1[j] = np.clip(off_1[j],self.bounds[j][0], self.bounds[j][1])
            # print("mutation",delta)
             
            off[suoyin, :] = off_1
            off[suoyin + 1, :] = off_2
            suoyin += 2
        for i in range(self.pop_size):
            for j in range(x_num):
                if self.flag[j]==1:
                    off[i][j]=round(off[i][j])
        chromo_offspring = off
        

        return chromo_offspring
    def main(self):

    #Main program starts here
        # pop_size = 20
        # max_gen = 150
        solution=self.first_population()
        # print(solution)
        # solution = [[random.uniform(low, high) for low, high in self.bounds] for i in range(0,self.pop_size)]
        for i in range(self.pop_size):
            for j in range(len(self.bounds)):
                if self.flag[j]==1:
                    solution[i][j]=round(solution[i][j])
        # print(solution)
        # solution=[[self.min_x[0]+(self.max_x[0]-self.min_x[0])*random.random() ,self.min_x[1]+(self.max_x[1]-self.min_x[1])*random.random()]for i in range(0,pop_size)]
        solution=np.array(solution)
        first_gen=solution 


        
         
        gen_no=0
        while(gen_no<self.max_gen):
            print('\n')
            print('gen_no:迭代次数',gen_no)
            function1_values = [self.obj1(solution[i])for i in range(0,self.pop_size)]
            
            function2_values = [self.obj2(solution[i])for i in range(0,self.pop_size)]
            
            CV_value=[self.constraint(solution[i])for i in range(0,self.pop_size)]
            # print('function1_values:',function1_values)
            # print('function2_values:', function2_values)
            non_dominated_sorted_solution = self.fast_non_dominated_sort(function1_values[:],function2_values[:],CV_value[:])
            # print('front',non_dominated_sorted_solution)
            # print("The best front for Generation number ",gen_no, " is")
            # for valuez in non_dominated_sorted_solution[0]:
            #     print("solution[valuez]",round(solution[valuez],3),end=" ")
        #     print("\n")
            crowding_distance_values=[]
            for i in range(0,len(non_dominated_sorted_solution)):
                crowding_distance_values.append(self.crowding_distance(function1_values[:],function2_values[:],non_dominated_sorted_solution[i][:]))
            # print("crowding_distance_values",crowding_distance_values)
            solution2 = solution[:]
             
            #Generating offsprings
            # while(len(solution2)!=2*pop_size):
            #     a1 = random.randint(0,pop_size-1)
            #     b1 = random.randint(0,pop_size-1)
            #     a2 = random.randint(0,pop_size-1)
            #     b2 = random.randint(0,pop_size-1)
            #     solution2=np.vstack((solution2,crossover(solution[a1],solution[b1])))
            solution2=np.vstack((solution,self.cross_mutation(solution,1,1,20,20)))
            # print(cross_mutation(solution,1,0.5,0.5,0.5))
            # print('solution2',solution2)
            function1_values2 = [self.obj1(solution2[i])for i in range(0,2*self.pop_size)]
            function2_values2 = [self.obj2(solution2[i])for i in range(0,2*self.pop_size)]
            CV_value=[self.constraint(solution2[i])for i in range(0,2*self.pop_size)]
            non_dominated_sorted_solution2 = self.fast_non_dominated_sort(function1_values2[:],function2_values2[:],CV_value[:])  #2*pop_size
            # print('non_dominated_sorted_solution2', non_dominated_sorted_solution2)
            # print("\n")
            crowding_distance_values2=[]
            for i in range(0,len(non_dominated_sorted_solution2)):
                crowding_distance_values2.append(self.crowding_distance(function1_values2[:],function2_values2[:],non_dominated_sorted_solution2[i][:]))
        
            # print('crowding_distance_values2',crowding_distance_values2)
            new_solution= []
            for i in range(0,len(non_dominated_sorted_solution2)):
            
                # non_dominated_sorted_solution2_1 = [index_of(non_dominated_sorted_solution2[i][j],non_dominated_sorted_solution2[i] ) for j in range(0,len(non_dominated_sorted_solution2[i]))]
                # print(non_dominated_sorted_solution2_1)
                # print(non_dominated_sorted_solution2)
                # print('non_dominated_sorted_solution2_1:',non_dominated_sorted_solution2_1)
                # front22 = sort_by_values(non_dominated_sorted_solution2_1[:], crowding_distance_values2[i][:])#对同层中的拥挤度进行排序
                front22 = self.sort_by_values(range(0,len(non_dominated_sorted_solution2[i])), crowding_distance_values2[i][:])#对同层中的拥挤度进行排序
                # print("front22",front22)

                front = [non_dominated_sorted_solution2[i][front22[j]] for j in range(0,len(non_dominated_sorted_solution2[i]))]
                # print(crowding_distance_values2[i])
                # print('front',front22)
                # print('front',front)
                front.reverse()
                for value in front:
                    new_solution.append(value)
                    if(len(new_solution)==self.pop_size):
                        break
                if (len(new_solution) == self.pop_size):
                    break
            solution = np.array([solution2[i] for i in new_solution])
            gen_no = gen_no + 1

            function1 = [i * -1 for  i in function1_values2]
            function2 = [j * -1  for j in function2_values2]
            feasible_index=[i for i in range(2*self.pop_size) if CV_value[i]==0]
            feasible_value1=[-1*function1_values2[i] for i in feasible_index]
            feasible_value2=[-1*function2_values2[i] for i in feasible_index]
            infeasible_index=[i for i in range(2*self.pop_size) if CV_value[i]>0]
            infeasible_value1=[-1*function1_values2[i] for i in infeasible_index]
            infeasible_value2=[-1*function2_values2[i] for i in infeasible_index]
            front_value1=[-1*function1_values2[i] for i in non_dominated_sorted_solution2[0]]
            front_value2=[-1*function2_values2[i] for i in non_dominated_sorted_solution2[0]]
            plt.scatter(feasible_value1,feasible_value2, marker=".",color=(0.2,0.5,gen_no/self.max_gen))
            # plt.scatter(infeasible_value1,infeasible_value2, marker="*",color='k')
            plt.scatter(front_value1,front_value2, marker="x",color=(0.2,0.5,gen_no/self.max_gen))
            
            # plt.show()
        # Lets plot the final front now
        # print(first_gen)
        # print(solution)
        print(CV_value)
        print(non_dominated_sorted_solution)
        # function1 = [i * -1  for i in function1_values]
        # function2 = [j * -1  for j in function2_values]
        function1 = [function1_values2[i] * -1  for i in new_solution]
        function2 = [function2_values2[i] * -1  for i in new_solution]
        # v1=1
        # v2=-1
        # index=0
        # for i in range(0,self.pop_size):
        #     if (function1[i]-v1)*(function2[i]-v2)>max_obj:
        #         max_obj=(function1[i]-v1)*(function2[i]-v2)
        #         index=i
        # if index==0 or index==self.pop_size:
        #     self.max_gen
        return [solution,function1,function2]

# def obj1(x):
#     value = -x[0]**2+2*x[1]**2
#     return value

# #Second function to optimize
# def obj2(x):
#     value = -(x[0]-2)**2-x[1]
#     return value

# def CV_value(x):
#     if x[0]<=55:
#         return 0
#     else :
#         return x[0]-3
    
    

# if __name__=='__main__':
#     NSGA=NSGA2(obj1,obj2,CV_value)
#     NSGA.main()