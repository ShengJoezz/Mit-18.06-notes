---
author: joe
rating: 8
content: 本文介绍了一种新的矩阵马尔可夫矩阵，其与差分方程的系数矩阵非常像，有特殊的性质与特殊的用处；同时介绍了标准正交基在傅里叶级数展开中的应用。
coverage: p24
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-07
---
# Section1.马尔科夫矩阵

## Subsection1.概念

马尔可夫矩阵是一种描述状态转移概率的矩阵，说白点就是基于上一个状态来预测下一个状态，例如两个地方每一年的人口转移；股票市场两只股票的关系等等。 ^2e6940

满足以下两个要求的矩阵即为马尔科夫矩阵（Markov matrix）：

1. 矩阵中的所有元素**大于等于**$0$；（因为马尔科夫矩阵与概率有关，而概率是非负的。）
2. 每一列的元素之和为$1$

## Subsection2.性质与证明

根据上述两个要求，可知马尔可夫矩阵会具有以下两个性质:

1. 马尔科夫矩阵必有特征值为$1$；
2. 其他的特征值的绝对值皆小于$1$。

### 证明:
#### 性质1:
##### 方法1:
因为马尔可夫矩阵每一列的元素之和为1.那么我们考虑构造如下情形：
$$x=\begin{bmatrix}1\\1\\\cdots\\1\end{bmatrix}\ \ \therefore x^TA=\begin{bmatrix}1&1&\cdots&1\end{bmatrix}=x^T$$
那么我们取转置会有:
$$A^Tx=x\therefore\lambda_{A^T}=1$$
因此我们现在需要证明$\lambda_{A^T}=\lambda_{A}$:
$$\begin{aligned}考虑我们对于A&的特征方程:det{(A-\lambda I)}=0\\而det{(A-\lambda I)}=0时，&由于我们知道转置的行列式等于转置前\\\therefore det(A-\lambda I)^T&=det(A^T-\lambda I)=0 \therefore \lambda_{A}=\lambda_{A^T}\end{aligned}$$
##### 方法2:
考虑$A-I$，此时每一列相加为0，那么我们直接把其余的所有行加到第一行，第一行就会变成零行，因此$det(A-I)=0$，即证。

#### 性质2:
涉及到更高阶的知识，我们可以证明马尔可夫矩阵的谱半径（即最大特征值的绝对值）不超过1。具体的证明方法有多种，这里不考虑证明了。

我们之前提到，[[#^2e6940|马尔可夫矩阵基于上一个状态来预测下一个状态]]，其实就是类似于[[Chapter20矩阵对角化，差分方程(A的乘方的简便计算)|差分方程的做法]]($u_{k+1}=Au_{k}$)。我们知道这种时候方程的解为:$u_k=A^ku_0=S\Lambda^kS^{-1}u_0=S\Lambda^kS^{-1}Sc=S\Lambda^kc=c_1\lambda_1^kx_1+c_2\lambda_2^kx_2+\cdots+c_n\lambda_n^kx_n$，由于我们提到的性质2:其他的特征值的绝对值皆小于$1$,那么于是在经过$k$次迭代，随着时间的推移，其他项都趋近于$0$，于是在$k\to\infty$时，有稳态$u_k=c_1x_1$，这也就是初始条件$u_0$的第$1$个分量。

## Subsection3.应用

接下来介绍马尔科夫矩阵的应用，我们用麻省和加州这两个州的人口迁移为例：

$\begin{bmatrix}u_{cal}\\u_{mass}\end{bmatrix}_{k+1}\begin{bmatrix}0.9&0.2\\0.1&0.8\end{bmatrix}\begin{bmatrix}u_{cal}\\u_{mass}\end{bmatrix}_k$，元素非负，列和为一。这个式子表示每年有$10\%$的人口从加州迁往麻省，同时有$20\%$的人口从麻省迁往加州。注意使用马尔科夫矩阵的前提条件是随着时间的推移，矩阵始终不变。

设初始情况$\begin{bmatrix}u_{cal}\\u_{mass}\end{bmatrix}_0=\begin{bmatrix}0\\1000\end{bmatrix}$，我们先来看第一次迁徙后人口的变化情况：$\begin{bmatrix}u_{cal}\\u_{mass}\end{bmatrix}_1=\begin{bmatrix}0.9&0.2\\0.1&0.8\end{bmatrix}\begin{bmatrix}0\\1000\end{bmatrix}=\begin{bmatrix}200\\800\end{bmatrix}$，随着时间的推移，会有越来越多的麻省人迁往加州，而同时又会有部分加州人迁往麻省。

计算特征值：我们知道马尔科夫矩阵的一个特征值为$\lambda_1=1$，则另一个特征值可以直接从迹算出$\lambda_2=0.7$。

计算特征向量：带入$\lambda_1=1$求$A-I$的零空间有$\begin{bmatrix}-0.1&0.2\\0.1&-0.2\end{bmatrix}$，则$x_1=\begin{bmatrix}2\\1\end{bmatrix}$，此时我们已经可以得出无穷步后稳态下的结果了。$u_{\infty}=c_1\begin{bmatrix}2\\1\end{bmatrix}$且人口总数始终为$1000$，则$c_1=\frac{1000}{3}$，稳态时$\begin{bmatrix}u_{cal}\\u_{mass}\end{bmatrix}_{\infty}=\begin{bmatrix}\frac{2000}{3}\\\frac{1000}{3}\end{bmatrix}$。注意到特征值为$1$的特征向量元素皆为正。

为了求每一步的结果，我们必须解出所有特征向量。带入$\lambda_2=0.7$求$A-0.7I$的零空间有$\begin{bmatrix}0.2&0.2\\0.1&0.1\end{bmatrix}$，则$x_2=\begin{bmatrix}-1\\1\end{bmatrix}$。

通过$u_0$解出$c_1, c_2$，$u_k=c_11^k\begin{bmatrix}2\\1\end{bmatrix}+c_20.7^k\begin{bmatrix}-1\\1\end{bmatrix}$，带入$k=0$得$u_0=\begin{bmatrix}0\\1000\end{bmatrix}=c_1\begin{bmatrix}2\\1\end{bmatrix}+c_2\begin{bmatrix}-1\\1\end{bmatrix}$，解出$c_1=\frac{1000}{3}, c_2=\frac{2000}{3}$。

另外，有时人们更喜欢用行向量，此时将要使用行向量乘以矩阵，其行向量各分量之和为$1$。
```ad-note
这里的难点我觉得是矩阵如何构造:
我们知道$A$是每一列的和是1，而第一列的元素是和后面的$x$向量的第一行相乘，因此也就是要说明矩阵的第$n$列是对$x$元素的第$n$行的划分，因此我们这里让第一列代表第一年的加州人数划分为第二年去加州与去麻州的人数。

```

# Section2.傅里叶级数

==我觉得非常精彩的一个part==

## Subsection1.重要的引入

我们在之前就已经讲过了任意一个向量$v$在标准正交基上的展开:
$$v=x_1q_1+x_2q_2+\cdots+x_nq_n$$
对于系数的求解是一个重点，这里我们巧妙地利用正交性与标准性，例如求解$x_1$，可以两边同乘以$q_1^T$，就会有:
$$q_1^Tv=x_1$$
我们用矩阵角度来考虑会有:
$$\Bigg[q_1\ q_2\ \cdots\ q_n\Bigg]\begin{bmatrix}x_1\\x_2\\\vdots\\x_n\end{bmatrix}=v\rightarrow Qx=v\rightarrow x=Q^{-1}v$$
这里注意到，又标准正交基构成的正交矩阵满足$Q^T=Q^{-1}$,因此原式改写为:
$$x=Q^Tv$$
此时对于$x$的每一个分量有:
$$x_i=q_i^Tv$$
## Subsection2.傅里叶级数

先写出傅里叶级数的展开式：

$$
f(x)=a_0+a_1\cos x+b_1\sin x+a_2\cos 2x+b_2\sin 2x+\cdots
$$

傅里叶发现，如同将向量$v$展开（投影）到向量空间的一组标准正交基中，在函数空间中，我们也可以做类似的展开。将函数$f(x)$投影在一系列相互正交的函数中。函数空间中的$f(x)$就是向量空间中的$v$；函数空间中的$1,\cos x,\sin x,\cos 2x,\sin 2x,\cdots$就是向量空间中的$q_1,q_2,\cdots,q_n$；不同的是，函数空间是无限维的而我们以前接触到的向量空间通常是有限维的。

但是这里我们要注意一点，我们选用的正交基是否正交呢？因此我们需要引入函数正交的概念:

我们知道正交的概念:
![[Chapter12正交向量与子空间#^4fe74c]]

因此我们要想说正交，就要考虑内积是什么:我们知道对于向量$v,w$的内积为$v^Tw=v_1w_1+v_2w_2+\cdots+v_nw_n=0$，也就是向量的每个分量之积再求和。而对于函数$f(x)\cdot g(x)$内积，同样的，我们需要计算两个函数的每个值之积而后求和，由于函数取值是连续的，所以函数内积为：
$$f^Tg=\int f(x)g(x)\mathrm{d}x$$
在本例中，由于傅里叶级数使用正余弦函数，它们的周期都可以算作$2\pi$，所以本例的函数点积可以写作$f^Tg=\int_0^{2\pi}f(x)g(x)\mathrm{d}x$。我来检验一个内积$\int_0^{2\pi}\sin{x}\cos{x}\mathrm{d}x=\left.\frac{1}{2}\sin^2x\right|_0^{2\pi}=0$。因此我们可以说明$sinx和cosx$是正交的。

```ad-note
title: 关于三角函数族的正交性
$\int_{-\pi}^{\pi}\cos(nx)\cos(mx)dx=\pi\delta_{mn},m,n\geq 1,$

$\int_{-\pi}^{\pi}\sin(nx)\sin(mx)dx=\pi\delta_{mn},m,n\geq 1$

$\int_{-\pi}^{\pi}\cos(nx)\sin(mx)dx=0;$
(这里的$\delta_{mn}$是克罗内克函数),而$\delta _{{ij}}=\left\{{\begin{matrix}1&(i=j)\\0&(i\neq j)\end{matrix}}\right.\,\!$

参见[傅里叶级数](https://zh.wikipedia.org/wiki/%E5%82%85%E9%87%8C%E5%8F%B6%E7%BA%A7%E6%95%B0)

```

==我们关注傅里叶级数的系数怎么求解==同样的，我们利用类似的手法:

$a_0$是容易得到的:$$\int_0^{2\pi}f(x)\mathrm{d}x =a_0\int_0^{2\pi} 1\times1\mathrm{d}x\rightarrow a_0=\frac{1}{2\pi}\int_0^{2\pi}f(x)\mathrm{d}x$$
对两边同时进行内积的求解:
$$\int_0^{2\pi}f(x)\cos x\mathrm{d}x=a_1\int_0^{2\pi}\cos^2x\mathrm{d}x=a_1\pi\rightarrow a_1=\frac{1}{\pi}\int_0^{2\pi}f(x)\cos x\mathrm{d}x$$

于是，我们把函数$f(x)$展开到了函数空间的一组标准正交基上。