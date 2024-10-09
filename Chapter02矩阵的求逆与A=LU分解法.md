---
author: joe
rating: 8
content: 基于行变换求解的Guass-Jordan方法求逆矩阵，同时介绍了A=LU分解方法
coverage: p3-p4
tags:
  - linear_algebra
  - math
  - class
---

# Section1. 矩阵何时会有逆

我们考虑一个矩阵:

$$A=
\begin{bmatrix}
1 & 2\\
3 & 6\\
\end{bmatrix}
$$
我们考虑是否会有$$A^{-1}A=I$$
事实上我们发现并没有，那我们从什么角度去考虑没有的? ^7403bf

* 后面的行列式会告诉我们，他的行列式是0
* 我们从[[Chapter01线性方程组的求解#Subsection5.矩阵乘法的几个视角|右乘是列变换]]的视角来看，注意到$A$ 中的两个列向量是共线的，所以我们无法通过他们的`linear combination`去实现形成一个$\begin{bmatrix}1\\0 \end{bmatrix}$ 的形式。 ^06ba1a

我们从第二点出发，其实也就是说明了一点，即:

$$\begin{bmatrix}1&2\\3&6\end{bmatrix}\cdot \begin{bmatrix}a\\b\end{bmatrix}=0\rightarrow AX=0$$
是**不应该有非零解**的。

那么我们怎么严谨地去说明呢？
$$if\ \ AX=0 \ \ \ \ \exists A^{-1}\ \ \ s.t.A^{-1}A=I$$
$$\therefore A^{-1}AX=0 \ \ X=0$$
可以发现，如果存在逆矩阵，则$X$ 只能是0向量.

# Section2. 对逆矩阵的求解方式(Gauss-Jordan方法)

由[[Chapter02矩阵的求逆与A=LU分解法#^06ba1a|求逆矩阵实际是列向量的线性组合]]这一角度来看，我们实际上是在做$AX=b\ \ b=\begin{bmatrix}1\\0\end{bmatrix}$ 的解方程过程，实际上我们在[[Chapter01线性方程组的求解#Subsecion4.利用初等矩阵进行的Gauss消元(行变换)]]已经提到过这个方程的解决方案，即通过$E_{12}$ 等初等矩阵进行行变换，使主元暴露。

但是，我们在进行逆矩阵求解的时候，逆矩阵的列很多，因此我们需要做多个方程的求解，那么我们可以把写在一起，即:

$$\begin{bmatrix}1&2|1&0\\
3&6|0&1\end{bmatrix}$$
进行消元的过程，直到把左侧行变换成为$I$ 的形式，此时右侧就是$A^{-1}$ .

$$\begin{bmatrix}1&0|1&3\\
0&1|2&7\end{bmatrix}$$

从理论上来分析，我们的左侧做了一个工作:

$$EA=I$$
那么:
$$E=A^{-1}$$
而右侧是:
$$IE=E=A^{-1}$$
忠实的记录了逆矩阵，因此可以得到$A^{-1}$ . ^5c9a4c

这就是Gauss-Jordan方法.

# Section3.转置的逆，逆的转置

$$if\ \ AA^{-1}=I$$
$$\because (AB)^T=B^TA^T$$
$$\therefore (A^{-1})^{T}A^T=I$$
那么说明:
$$(A^T)^{-1}=(A^{-1})^T$$
这一点就说明转置与逆运算是可交换的运算.

# Section4.A=LU的分解

首先解释一下$A=LU$ 是什么，$U$ 在[[Chapter01线性方程组的求解#^5816e6|行变换消元的最后结果]]中已经有提到，实际上是`Upper triangular`(上三角矩阵)，而$L$ 是`Lower triangular`(下三角矩阵).

我们举几个例子来引入:
## Subsection1.二阶矩阵的简单例子

首先我们考虑某个简单的二阶矩阵，例如:

$$A=\begin{bmatrix}2 &1\\8&7\end{bmatrix}$$
如果我们要把他消成$U$ ，显然我们利用$E_{21}$ :

$$E_{21}=\begin{bmatrix}1&0\\-4&1\\\end{bmatrix}$$
如此的话，我们可以调整成为:
$$E_{21}A=U=\begin{bmatrix}2&1\\0&3\end{bmatrix}$$
那么这里我们只需要考虑取${E_{21}}^{-1}=L$ ,就可以得到$L$ 了:

$$L=\begin{bmatrix}1&0\\4&1\end{bmatrix}$$
这是个容易的例子，我们在求解$E_{21}$ 和求解$L$ 的难度几乎相差不大.

## Subsection2.三阶矩阵陡增的难度

这时，如果我们考虑一个三阶矩阵的例子，我们知道，利用[[Chapter01线性方程组的求解#Subsecion4.利用初等矩阵进行的Gauss消元(行变换)|Guass消元法]]，我们可以顺利地得到:

$$E_{32}E_{31}E_{21}A=U$$
那么我们知道我们可以利用求逆，得到:

$$L={E_{21}}^{-1}{E_{31}}^{-1}{E_{32}}^{-1}$$
这里其实我们就会发现，同样是乘法，为什么我们需要用${E_{21}}^{-1}{E_{31}}^{-1}{E_{32}}^{-1}$  而非$E_{32}E_{31}E_{21}$ 的形式，例如我们考虑两个如下的初等矩阵:

$$E_{32}\times E_{21}=\begin{bmatrix}1&0&0\\0&1&0\\0&-5&1\end{bmatrix}\times\begin{bmatrix}1&0&0\\-2&1&0\\0&0&1\\\end{bmatrix}=\begin{bmatrix}1&0&0\\-2&1&0\\10&-5&1\end{bmatrix}$$
我们会发现，这里竟然出现了$10$ 这个并不属于我们初等变化里面的数字!

事实上，这里就是因为我们在$E_{21}$ 中第二行有变形，而$E_{32}$ 又要用到第二行，导致变形变得更大了。

那假如我们考虑他们的逆矩阵,有:

$${E_{21}}^{-1}{E_{32}}^{-1}=\begin{bmatrix}1&0&0\\2&1&0\\0&0&1\\\end{bmatrix}\times\begin{bmatrix}1&0&0\\0&1&0\\0&5&1\end{bmatrix}=\begin{bmatrix}1&0&0\\2&1&0\\0&5&1\end{bmatrix}$$
我们会发现这里忠实地保留了初等矩阵中消元乘数的值，而没有去增加新的不必要的内容.

事实上，我们最后得到的$L=\begin{bmatrix}1&0&0\\-i&1&0\\-j&-k&1\end{bmatrix}$ 正是消元乘数的值的负数，并且求解也不困难.

## Subsection3.我们为什么要用A=LU？

就[[Chapter02矩阵的求逆与A=LU分解法#Subsection2.三阶矩阵陡增的难度]]的末尾我们可以知道，$L$ 容易得到，并且不会破坏消元乘数(事实上这两者可能是一回事，毕竟就是因为不会破坏我们才容易得到),我们只需要关注$U$ ，就比较容易得到$L$ 
