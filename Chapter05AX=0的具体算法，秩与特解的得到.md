---
author: joe
rating: 8
content: Ax=0求解的具体算法，特解怎么书写，为什么这样书写，列秩=行秩=矩阵秩是为什么
coverage: p7
tags:
  - linear_algebra
  - math
  - class
---
# Section1.概念引入

我们之前对于$Ax=0$的计算只考虑了宏观层面的，对于具体的算法我们这个part来进行介绍，首先我们考虑这个矩阵:

$$A=\begin{bmatrix}1&2&2&2\\2&4&6&8\\3&6&8&10\end{bmatrix}$$
我们采用[[Chapter01线性方程组的求解#Subsecion3.消元法]]把$A$转化为$U$，即:

$$U=\begin{bmatrix}1&2&2&2\\0&0&2&4\\0&0&0&0\end{bmatrix}$$
那么这里我们注意到几个概念:

* 主元:这里我们这个矩阵有两个主元，分别是$a_{11}$与$a_{23}$即阶梯型矩阵的阶梯处.
* 主列与自由列:几列对应几个变量，有主元的列成为主列，其余的列即自由列，都是自由变量.
* 秩:主元的个数.

# Section2.一些新的视角

## Subsection1.关于$Ax=0$与$Ux=0$

我们在[[#Section1.概念引入]]中其实已经在把求解$Ax=0$转换为求解$Ux=0$了，那么这两个等式等价吗？

结论是显然肯定的，因为我们做行变换并不破坏这些等式本身，这些行变换前后本身就是等价的，换句话来说，这两个等式:

**零空间相同，行变换不会改变零空间**

## Subsection2.关于行秩，列秩，以及秩

在这个例子里，秩=主元个数=2，且我们注意到列向量与行向量的线性无关的向量数目均是3，其实我们就可以还注意到:

**列秩=行秩=矩阵秩** ^7af8fa

## Subsection3.转置的矩阵的秩与原矩阵的关系

r($A$)=r($A^T$)

矩阵的主元的个数与其转置的主元的个数相同。

主元的概念是伴随消元而生的。假设矩阵 $A$ 消元后所得的主元个数为 $a$，那么就表示矩阵 $A$ 有 $a$ 行线性无关，有 $a$ 列线性无关；假设矩阵 $A^T$ 消元后所得的主元个数为 $b$，那么就表示矩阵 $A^T$ 有 $b$ 行线性无关，有 $b$ 列线性无关。

而由矩阵的转置的定义，我们知道，矩阵 $A$ 的行就是矩阵 $A^T$ 的列，因此 $a=b$。

换句话说，$A$ 的主元个数 = $A$ 矩阵线性无关的列的个数 = $A^T$ 矩阵线性无关的行的个数 = $A^T$ 的主元个数。
### tip:做行变换为什么可以知道列秩?

行秩等于矩阵秩我想是比较好理解的，因为我们是在对行向量进行变换当然会发现其中的线性相关性，那么我们做行变换为什么会反应列秩呢？这里我们做一个简单的证明:

$$A=\begin{bmatrix}a&c&k_1a+k_2c\\b&d&k_1b+k_2d\\e&f&k_1e+k_2f\end{bmatrix}$$
我们经过行变换可以得到:

$$A=\begin{bmatrix}a&c&k_1a+k_2c\\0&\frac{ad-bc}{a}&k_2\frac{ad-bc}{a}\\0&\frac{af-ec}{a}&k_2\frac{af-ec}{a}\end{bmatrix}$$
这个时候我们注意到，我们做了行变换，但是我们的列三($\beta_3$)仍旧满足$\beta_3=k_1\beta_1+k_2\beta_2$，说明:

* 行变换并没有改变列的线性组合
* 即便我们不能确定行向量的线性关系是什么样子的，但是列向量的秩也能说明行向量的秩

# Section3.特解与方程的解

我们现在已经得到了
$$U=\begin{bmatrix}1&2&2&2\\0&0&2&4\\0&0&0&0\end{bmatrix}$$
我们可以做一些进一步的化简，把主元所在的列全部清零(除开主元处)

$$R=\begin{bmatrix}1&2&0&-2\\0&0&1&2\\0&0&0&0\end{bmatrix}$$
这个$R$我们又称为$RREF$(Reduced Row Echelon Form),即约化行阶梯形式，在$matlab$中也有使用.

## Subsection1.方程解向量的数量

我们这里主元是$x_1,x_3$，自由变量是$x_2,x_4$,那么我们的解向量可以写成如下形式:

$$X=\begin{bmatrix}2x_4-2x_2\\x_2\\-2x_4\\x_4\end{bmatrix}=x_2\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}+x_4\begin{bmatrix}2\\0\\-2\\1\end{bmatrix}$$
显然，我们解向量的数量依赖于自由变量的数量，而**自由变量的数量是$n(变量数目)-r(秩)$.**

## Subsection2.方程的特解

我们在[[#Subsection1.方程解向量的数量]]中已经写出了两个解向量$\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}$与$\begin{bmatrix}2\\0\\-2\\1\end{bmatrix}$,这就是特解，当然我们还可以用另一个角度来求解特解，即:
令自由变量$\begin{bmatrix}x_2\\x_4\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix}或\begin{bmatrix}0\\1\\\end{bmatrix}$这样得到两组解，分别为特解1与特解2.

这种做法是有讲究的，我们令其中一个自由变量为0，另一个为1，是为了专注于观察某一个特定的变量对方程组解的影响. ^1e2bd1

# Section3.对于$RREF$矩阵和解的理解

这是我们得到的$RREF$矩阵:
$$R=\begin{bmatrix}1&2&0&-2\\0&0&1&2\\0&0&0&0\end{bmatrix}$$
我们把自由列$x_2$和主列$x_3$交换一下便于观察:
$$R=\begin{bmatrix}1&0&2&-2\\0&1&0&2\\0&0&0&0\end{bmatrix}$$
我们的解向量放到一个矩阵里并且同样交换$x_2$和$x_3$的位置，得到:
$$\begin{bmatrix}-2&2\\0&-2\\1&0\\0&1\end{bmatrix}$$
我们发现，主元所在的列形成了一个$I$放到了最下，而自由列则是取了相反数放到了上面. ^bbe1c6

这也是很多考研的书上告诉我们的求解线性方程组的方式.

## Prove:为什么可以这样求解线性方程组？

我们注意到，我们始终可以把$RREF$矩阵化简为如下形式:

$$R=\begin{bmatrix}I&F\\0&0\end{bmatrix}\ \ F指的是自由列形成的矩阵$$
因此我们考虑的:
$$RX=0$$
可以知道:
$$X=\begin{bmatrix}-F\\I\end{bmatrix}$$
那么这里我们可以对为什么$X=\begin{bmatrix}-F\\I\end{bmatrix}$进行多一点的说明:

我们考虑一个零空间矩阵，$$X=\begin{bmatrix}x_1&x_2&...&x_{n-r}\end{bmatrix}\ \ 其中x_k是RX=O的特解$$
我们可以把$X$写为:

$$X=\begin{bmatrix}X_{pivot}\\X_{free}\end{bmatrix}$$
那么我们按照[[#^1e2bd1|令自由变量一个为1其余都为0]]的手法，相当于让$X_{free}=I$(真是一个让人惊喜的巧合！)，那么我们顺利地得到:

$$X=\begin{bmatrix}-F\\I\end{bmatrix}$$
## 势利一点的结论:方程化简到了$RREF$的形式，怎么得到最后的特解？

$e.g.$  

$$R=\begin{bmatrix}1&a&0&c\\0&b&1&d\\0&0&0&0\end{bmatrix}$$
step1:自由变量是$x_2,x_4$，那么先把$x_2,x_4$的位置占住，用$\begin{bmatrix}x_2\\x_4\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix}或\begin{bmatrix}0\\1\\\end{bmatrix}$
step2:占住之后，最前的一列自由列取相反数与$\begin{bmatrix}1\\0\end{bmatrix}$组合，后面的相应组合.




