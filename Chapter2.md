# 高等数值计算第二章复习
每一章，包括一些公式都给整理下来便于记忆。复习完这一遍之后就可以去整理一个公式纸了。

## 知识结构图

![image-20231227142331672](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-1)

## 拉格朗日插值多项式

已知$y=f(x)$在区间$[x_0,x_n]$上有定义及在$n+1$个节点$x_0<x_1<...<x_n$的函数值$y_j=f(x_j)(j=0,1,...,n)$

要求$n$次插值多项式$L_n(x)$，使他满足
$$
L_n(x_j)=y_j,(j=0,1,...,n)
$$


### 基函数法求解

$L_n(x)$表示为已知节点函数值的基函数组合形式：
$$
L_n(x)=\sum _{k=0}^{n} y_{k} l_{k}( x)
$$
其中，组合系数被称为$y_k$，而$l_k(x)$则称为$n$次插值基函数，满足下面条件：

1. $l_k(x)(k=0,1,...,n)$是不超过$n$的多项式函数

2. 在节点$x_k(k=0,1,...,n)$处满足：
   $$
   l_{k}( x) =\begin{cases}
   1, & k=j\\
   0, & k\neq j
   \end{cases}( j,k=0,1,...,n)
   $$

推导过程在此处，核心思路就是做一个待定系数法，把那个系数给求出来：

![image-20231227144728461](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-2)

所以n次拉格朗日插值多项式$L_n(x)$为：
$$
L_n(x)=\sum _{k=0}^{n} y_{k} \frac{\omega_{n+1}(x)}{(x-x_k)\omega'_{n+1}(x_k)}
$$

### 插值余项

若在$[a,b]$上用$L_n(x)$近似$f(x)$，则其截断误差为$R=f(x)-L_n(x)$，也称为插值多项式的余项，记作$R_n(x)$

设$f^{(n)}(x)$在$[a,b]$上连续，$f^{(n+1)}(x)$在$[a,b]$上存在，节点$a\le x_0<x_1<...<x_n\le b$，则$L_n(x)$是满足拉格朗日插值条件的多项式，则对任何$x\in [a,b]$，有：
$$
R_n(x)=f(x)-L_n(x)=\frac{f^{n+1}(\xi)}{(n+1)!}\omega_{n+1}(x)
$$
这里$\xi \in [a,b]$且依赖于$x$

推导如下图：

![image-20231227150158021](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-3)

### 罗尔定理以及其Generalized形式

![image-20231227150248379](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-4)

### 余项的误差限

余项的误差限为：$R_n(x)\le \frac{M_{n+1}}{(n+1)!}|\omega_{n+1}(x)|$

### 例题

例1：

![image-20231227151748498](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-5)

例1就是按照拉格朗日插值的方法求出求出来各个$l_n(x)$

例2：

![image-20231227151922235](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-6)

![image-20231227152147474](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-7)

这个例题就是直接用二次插值的多项式来直接解

例3

![image-20231227152212869](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-8)

![image-20231227152232996](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-9)



因为要求的是这个x，而不是函数值，所以这里就做一个从函数值到自变量的一个反函数，并且利用这个反函数做插值。

例4

![image-20231227152729218](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-10)

第二题就是利用了第一题的性质去做的。

例5

![image-20231227153035064](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-11)

利用了之前的误差余项性质。



## 均差 & 牛顿插值多项式

### 均差定义与性质

均差的定义：

- 称$f[x_0,x_k]=\frac{f(x_k)-f(x_0)}{x_k-x_0}$为函数$f(x)$关于点$x_0,x_k$的一阶均差
- 称$f[x_0,x_1,x_k]=\frac{f[x_0,x_k]-f[x_0,x_1]}{x_k-x_1}$为函数$f(x)$关于点$x_0,x_1,x_k$的二阶均差
- 称$f[x_0,x_1,...,x_k]=\frac{f[x_0,...,x_{k-2},x_k]-f[x_0,x_1,...,x_{k-2},x_{k-1}]}{x_k-x_{k-1}}$为函数$f(x)$关于点$x_0,x_1,...,x_k$的$k$阶均差

均差的性质：

- $k$阶均差可表示为函数值$f(x_0),...,f(x_k)$的线性组合
  $$
  f[x_0,x_1,...,x_k]=\sum _{j=0}^{k} \frac{f(x_j)}{(x_j-x_0)...(x_j-x_{j-1})(x_j-x_{j+1})...(x_j-x_{k})}
  $$

- $f[x_0,x_1,...,x_k]=\frac{f[x_1,...,x_k]-f[x_0,...,x_{k-1}]}{x_k-x_{0}}$

- 若$f(x)$在$[a,b]$上存在$n$阶乘导数，且节点$x_0,..,x_n\in[a,b]$，则有：
  $$
  f[x_0,...,x_n]=\frac{f^{(n)}(\xi)}{n!}, \xi\in[a,b]
  $$
  

### 均差计算表

![image-20231227154643151](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-12)

![image-20231227154711790](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-13)

### 牛顿插值公式

通过均差的定义，用迭代的方法可以求得牛顿插值公式：
$$
N_n(x)=f(x_0)+f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+...+f[x_0,x_1,...,x_n](x-x_0)...(x-x_{n-1})
$$
而其余项则为：
$$
R_n(x)=f[x,x_0,...,x_n](x-x_0)...(x-x_n)
$$
例1

![image-20231227155055714](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-14)

这里就是先求出均差表，然后逐项把牛顿插值公式给写出来。

例2

![image-20231227155202568](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-15)

![image-20231227155231534](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-16)

这里就是直接用均差表求出来牛顿公式，然后再根据误差公式来算出他的误差



## 埃尔米特插值多项式

### 定义

定义：在节点上的函数值相等，并且对应的导数值也相等，甚至要求高阶要求也相等。

均差的性质：

![image-20231227155639892](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-17)

拓展到n次之后，则有：

![image-20231227155710697](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-18)

### 例题

![image-20231227155754429](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-19)

![image-20231227155818965](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-20)

![image-20231227155912067](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-21)

这里的第一种解决方法就是回到了定义，用定义的方式去解决问题。而第二种解决方法则是通过待定系数来解决问题。



例2

![image-20231227160016880](C:\Users\84865\Desktop\mypro\Advanced-Numerical-Calculation-Methods-Review\pics\2-22)

![image-20231227160032068](C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227160032068.png)

同样，第一个方法是用定义去算，第二个方法则是用待定系数法去算。



## 其它插值公式

### 分段线性插值

![image-20231227160439192](C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227160439192.png)

它的插值余项则有：

![image-20231227161144420](C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227161144420.png)

### 分段三次Hermite插值

![image-20231227161232933](C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227161232933.png)

### 三次样条插值

<img src="C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227161435315.png" alt="image-20231227161435315" style="zoom:150%;" />

![image-20231227161506414](C:\Users\84865\AppData\Roaming\Typora\typora-user-images\image-20231227161506414.png)