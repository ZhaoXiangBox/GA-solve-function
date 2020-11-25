/*******************************************************************************************    
 * 
 * 遗传算法求解步骤
 * 1.使用随机方法产生一个有N个染色体的初始种群
 * 2.对群体中的每一个染色体，计算其适应值
 * 3.从群体中随机选择一些染色体构成一个新的群体（常用方法:轮盘赌选择，锦标赛选择）
 * 4.以一定概率进行交叉(单点交叉、多点交叉)
 * 5.以一定概率进行变异
 * 6.返回第2步，以一定的迭代次数进行循环
 * 
 * 题目：求解一元函数  f(x) = x * sin(10 * PI * x) + 2.0   在[-1,2]中的最大值
 * Author： Lee  2020.11.25.
 *******************************************************************************************/
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>
using namespace std;

#define  PI  3.14159265

const double pcross = 0.7;												  //交叉概率
const double pmutate = 0.001;										//变异概率

// 取值范围的区间长度为3，pow(2,21) < 3 *1000000<pow(2,22)，因此编码的二进制最少为22位，即22个0和1
const int Gene_length = 22;											    //22位的染色体
const int Generations = 500;									  		//更迭代数
const int Size = 500;									     					 //种群规模
double bestval;												   						 //适应值最大值

typedef struct node {									//染色体结构体
	bool Chromosome[Gene_length];
}node;

node bestChromosome;                                       //记录最优个体
node group[Size];									     //记录种群中的个体的数组
node temp[Size];					  					 //记录种群中的个体的临时数组

/******  对单个染色体随机编码赋值 ******/
void Code_Chromosome(node& c) 
{								
	for (int i = 0; i < Gene_length; i++)
	{
		c.Chromosome[i] = rand() % 2;
	}
}

/*****  二进制解码操作  *****/
void decode(node& c, double& x) 
{		
	double num = 4194394;						   //即2的22次方
	double tem = 0;
	for (int i = 0; i < Gene_length;i++)
	{
		tem += c.Chromosome[i] * pow(2, i);		
	}
	x = (3 / num * tem)-1;
}

/*****  目标函数  *****/
double f(double x) 
{																					
	return x * sin(10 * PI * x) + 2.0;
}

/*****  适应度函数   （解码，换算成十进制求数值，并返回数值）*****/
double fitness(node& c) 
{								
	double x;
	decode(c, x);
	return f(x);
}

/*****  交叉操作  ( 实际就是两个数组元素互换的实现) *****/
void cross(node& c1, node& c2, int point) 
{	
	node c3 = c1;      // 元素互换，依赖于有一个中间元素暂时存储其中的一个元素
	for (int i = 0; i < Gene_length - point; i++)
	{
		c1.Chromosome[point + i ] = c2.Chromosome[point + i ];
	}
	for (int j = 0; j < Gene_length - point; j++) 
	{
		c2.Chromosome[point + j ] = c3.Chromosome[point + j ];
	}
}

/***** 变异操作 *****/
void mutate(node& c) 
{
	int i = rand() % Gene_length;        // rand()函数会输出一个伪随机数，这里对22求余数后，即得出被变异的某个 基因角标。
	c.Chromosome[i] = !c.Chromosome[i];
}

/*****  产生0到1的随机小数  *****/
double inline rand0()
{
	return rand() % 10000 / 10000.0;
}

/*****  选择操作  *****/
void select(node group[Size]) 				
{
	double fitnessval[Size];
	double sum = 0;
	double avgfitness[Size];
	int id[Size];
	for (int i = 0; i < Size; i++) 
	{
		fitnessval[i] = fitness(group[i]);
	}
	for (int i = 0; i < Size; i++) 						//适应度总和
	{
		sum += fitnessval[i];
	}

	for (int i = 0; i < Size; i++) 
	{
		avgfitness[i] = fitnessval[i] / sum;
	}

	for (int i = 1; i < Size; i++)  					//适应度累加
	{
		avgfitness[i] += avgfitness[i - 1];
	}

	for (int i = 0; i < Size; i++)              		 //轮盘赌选择法
	{
		double rannum = rand0();				//自定义函数，产生0到1随机数
		int j;
		for (j = 0; j < Size - 1; j++) 
		{
			if (rannum < avgfitness[j]) 
			{
				id[i] = j;
				break;
			}
		}
		if (j == Size - 1)
		{
			id[i] = j;
		}
	}

	for (int i = 0; i < Size; i++) //将新个体替换旧个体
	{
		temp[i] = group[i];
	}

	for (int i = 0; i < Size; i++)
	{
		group[i] = temp[id[i]];
	}
}

/***** 取得最优个体对应的位置  *****/
int GetBest(node group[Size], double& x, double& number) 
{
	double fitnessval[Size];
	for (int i = 0; i < Size; i++)
	{
		fitnessval[i] = fitness(group[i]);
	}
	int id = 0;
	for (int i = 1; i < Size; i++) 
	{
		if (fitnessval[i] > fitnessval[id]) 
		{
			id = i;
		}
	}
	decode(group[id], x);
	number = f(x);
	return id;
}

/*****  遗传算法流程 选择 交叉 变异 求当代最优*****/
void GA_Compute(double& x, double& number) 
{
	for (int i = 0; i < Size; i++) 
	{
		Code_Chromosome(group[i]);
	}
	bestChromosome = group[GetBest(group, x, bestval)];

	for (int i = 0; i < Generations; i++)         // 执行 Generations 代
	{
		select(group);											// 选择操作，
		int p = rand() % Gene_length;
		for (int j = 0, pre = -1; j < Size; j++)   // 根据概率交叉		
		{
			if (rand0() < pcross)
			{
				if (pre == -1)
					pre = j;
				else
				{
					cross(group[pre], group[j], p);
					pre = -1;
				}
			}
		}
		for (int k = 0, pre = -1; k < Size; k++)  //根据概率进行变异
		{
			if ((rand0() < pmutate))
			{
				mutate(group[k]);
			}
		}
		GetBest(group, x, number);
		cout<<"第"<<i<<"代中"<<"最优X值为:"<<x<<"\t对应的函数值为:"<<f(x)<<endl; 
	}

}

int main()
{
	srand((unsigned)time(0));      //产生随机数种子
	double x;
	double max;
	GA_Compute(x, max);
	return 0;
}