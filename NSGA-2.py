# Program Name: NSGA-II.py
# author；lxy
# Time:2023/03/11
 
#Importing required modules
import math
import random
import matplotlib.pyplot as plt
import numpy as np
import sys
#First function to optimize

class nsga2():
    def __init__(self):
    #Initialization
        self.min_x=[-55,-55]
        self. max_x=[55,55]
    def obj1(self,x):
        value = -x[0]**2
        return value
    
    #Second function to optimize
    def obj2(self,x):
        value = -(x[0]-2)**2
        return value
    
    #Function to find index of list,且是找到的第一个索引
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
            values[self.index_of(min(values),values)] = math.inf
        return sorted_list
    
    #Function to carry out NSGA-II's fast non dominated sort
    def fast_non_dominated_sort(self,values1, values2):
        S=[[] for i in range(0,len(values1))]   #len(values1)个空列表
        front = [[]]
        n=[0 for i in range(0,len(values1))]
        rank = [0 for i in range(0, len(values1))]
        #将front0全部整理出来了，并未对front1-n等进行整理
        for p in range(0,len(values1)):
            S[p]=[]
            n[p]=0
            for q in range(0, len(values1)):
                if (values1[p] > values1[q] and values2[p] > values2[q]) or (values1[p] >= values1[q] and values2[p] > values2[q]) or (values1[p] > values1[q] and values2[p] >= values2[q]):
                    if q not in S[p]:
                        S[p].append(q)
                elif (values1[q] > values1[p] and values2[q] > values2[p]) or (values1[q] >= values1[p] and values2[q] > values2[p]) or (values1[q] > values1[p] and values2[q] >= values2[p]):
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
            for k in range(2,lenth-1):
                distance[k] = distance[k]+ (values1[sorted1[k+1]] - values1[sorted1[k-1]])/(max(values1)-min(values1))
                # print("/n")
                print("k:",k)
                print("distance[{}]".format(k),distance[k])
            for k in range(2,lenth-1):
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
                    off_1[j] = np.clip(off_1[j], self.min_x[j], self.max_x[j])
                    off_2[j] = np.clip(off_2[j], self.min_x[j], self.max_x[j])

        

            # Polynomial mutation
            if random.random() < pm:
                u2 = np.random.rand(x_num)
                delta = np.zeros(x_num)
                for j in range(x_num):
                    if u2[j] < 0.5:
                        delta[j] = (2 * u2[j]) ** (1 / (yita2 + 1)) - 1
                    else:
                        delta[j] = 1 - (2 * (1 - u2[j])) ** (1 / (yita2 + 1))

                    off_1[j] += delta[j]
                    off_1[j] = np.clip(off_1[j], self.min_x[j],self.max_x[j])

        
            off[suoyin, :] = off_1
            off[suoyin + 1, :] = off_2
            suoyin += 2

        chromo_offspring = off
        return chromo_offspring
    def main(self):

    #Main program starts here
        pop_size = 16
        max_gen = 150
        

        
        solution=[[self.min_x[0]+(self.max_x[0]-self.min_x[0])*random.random() ,self.min_x[1]+(self.max_x[1]-self.min_x[1])*random.random()]for i in range(0,pop_size)]
        solution=np.array(solution)
        print('solution',solution.shape)


        

        gen_no=0
        while(gen_no<max_gen):
            print('\n')
            print('gen_no:迭代次数',gen_no)
            function1_values = [self.obj1(solution[i])for i in range(0,pop_size)]
            
            function2_values = [self.obj2(solution[i])for i in range(0,pop_size)]
            print('function1_values:',function1_values)
            print('function2_values:', function2_values)
            non_dominated_sorted_solution = self.fast_non_dominated_sort(function1_values[:],function2_values[:])
            print('front',non_dominated_sorted_solution)
            # print("The best front for Generation number ",gen_no, " is")
            # for valuez in non_dominated_sorted_solution[0]:
            #     print("solution[valuez]",round(solution[valuez],3),end=" ")
        #     print("\n")
            crowding_distance_values=[]
            for i in range(0,len(non_dominated_sorted_solution)):
                crowding_distance_values.append(self.crowding_distance(function1_values[:],function2_values[:],non_dominated_sorted_solution[i][:]))
            print("crowding_distance_values",crowding_distance_values)
            solution2 = solution[:]
            print(solution.shape)
            #Generating offsprings
            # while(len(solution2)!=2*pop_size):
            #     a1 = random.randint(0,pop_size-1)
            #     b1 = random.randint(0,pop_size-1)
            #     a2 = random.randint(0,pop_size-1)
            #     b2 = random.randint(0,pop_size-1)
            #     solution2=np.vstack((solution2,crossover(solution[a1],solution[b1])))
            solution2=np.vstack((solution,self.cross_mutation(solution,1,0.5,2,2)))
            # print(cross_mutation(solution,1,0.5,0.5,0.5))
            print('solution2',solution2)
            function1_values2 = [self.obj1(solution2[i])for i in range(0,2*pop_size)]
            function2_values2 = [self.obj2(solution2[i])for i in range(0,2*pop_size)]
            non_dominated_sorted_solution2 = self.fast_non_dominated_sort(function1_values2[:],function2_values2[:])  #2*pop_size
            print('non_dominated_sorted_solution2', non_dominated_sorted_solution2)
            # print("\n")
            crowding_distance_values2=[]
            for i in range(0,len(non_dominated_sorted_solution2)):
                crowding_distance_values2.append(self.crowding_distance(function1_values2[:],function2_values2[:],non_dominated_sorted_solution2[i][:]))
        
            print('crowding_distance_values2',crowding_distance_values2)
            new_solution= []
            for i in range(0,len(non_dominated_sorted_solution2)):
            
                # non_dominated_sorted_solution2_1 = [index_of(non_dominated_sorted_solution2[i][j],non_dominated_sorted_solution2[i] ) for j in range(0,len(non_dominated_sorted_solution2[i]))]
                # print(non_dominated_sorted_solution2_1)
                # print(non_dominated_sorted_solution2)
                # print('non_dominated_sorted_solution2_1:',non_dominated_sorted_solution2_1)
                # front22 = sort_by_values(non_dominated_sorted_solution2_1[:], crowding_distance_values2[i][:])#对同层中的拥挤度进行排序
                front22 = self.sort_by_values(range(0,len(non_dominated_sorted_solution2[i])), crowding_distance_values2[i][:])#对同层中的拥挤度进行排序
                print("front22",front22)
                front = [non_dominated_sorted_solution2[i][front22[j]] for j in range(0,len(non_dominated_sorted_solution2[i]))]
                print('front',front)
                front.reverse()
                for value in front:
                    new_solution.append(value)
                    if(len(new_solution)==pop_size):
                        break
                if (len(new_solution) == pop_size):
                    break
            solution = np.array([solution2[i] for i in new_solution])
            gen_no = gen_no + 1
        #
        # Lets plot the final front now
        function1 = [i * -1 for i in function1_values]
        function2 = [j * -1 for j in function2_values]
        plt.xlabel('Function 1', fontsize=15)
        plt.ylabel('Function 2', fontsize=15)
        plt.scatter(function1, function2)
        plt.show()

if __name__=='__main__':
    NSGA=nsga2()
    NSGA.main()