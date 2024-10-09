---
author: joe
rating: 8
content: 本文介绍了能够对角化的矩阵，介绍了他的定义，使用条件，证明与使用场景，同时说明了差分方程中的应用，并介绍了其在实际问题场景中的应用
coverage: p22
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-06
---
# Section1.矩阵对角化

上一章我们讲述了特征值$\lambda$与特征向量$x$的求法，那这里我们可以考虑一些特殊的方法。

我们在之前介绍了$A$的很多拆解方法，例如$A=LU$(利用[[Chapter02矩阵的求逆与A=LU分解法|行变换矩阵的求逆]]来拆解的)与$A=QR$(利用[[Chapter15正交矩阵与Gram-Schmidt正交化|施密特正交化]]来拆解)的方式，这里我们介绍一个新的拆解方式，即:$$A=S\Lambda S^{-1}$$
## Subsection1.证明

可以知道$S=\Bigg[x_1x_2\cdots x_n\Bigg]$,而$\lambda=\begin{bmatrix}\lambda_1&0&\cdots&0\\0&\lambda_2&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&\lambda_n\end{bmatrix}$,那么:
$$\begin{aligned}AS&=A\Bigg[x_1x_2\cdots x_n\Bigg]=\Bigg[\lambda_1x_1\ \ \ \lambda_2x_2\cdots \ \ \ \lambda_nx_n\Bigg]\\&=\Bigg[x_1x_2\cdots x_n\Bigg]\begin{bmatrix}\lambda_1&0&\cdots&0\\0&\lambda_2&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&\lambda_n\end{bmatrix}=S\lambda\\ &AS=S\lambda\rightarrow S^{-1}AS=\lambda\ or\ A=S\Lambda S^{-1}\end{aligned}$$
注意$S$的位置，只要记得推导过程先求了$AS$，那对角化之后右侧肯定是$S^{-1}$.

## Subsection2.使用的条件

我们根据式子形式很明显就可以看出来，由于有$S^{-1}$的出现，肯定要求$S$是满秩的，因此特征向量是必须互相线性无关的，因此有两种情况:

* 如果一个矩阵有$n$个互不相同的特征值（即没有重复的特征值），则该矩阵具有$n$个线性无关的特征向量，因此该矩阵可对角化。 ^fb51de
* 如果一个矩阵的特征值存在重复值，则该矩阵可能具有$n$个线性无关的特征向量。比如取$10$阶单位矩阵，$I_{10}$具有$10$个相同的特征值$1$,但是会有十个线性无关的特征向量，那也可以对角化，但是如果是[[Chapter19特征值和特征向量#^a8c6d5|退化矩阵]]，某个$n$重特征值只有小于$n$个线性无关的特征向量，那就不可以对角化。 ^45af7c

## Subsection3.应用

这个在关于$A^n$的计算与性质探究中很有用:
$$\begin{aligned}A^2=S\Lambda S^{-1}&S\Lambda S^{-1}=S\Lambda^2S^{-1}\\A^n=&S\lambda^{n}S^{-1}\end{aligned}$$
因此在求解$A^n$的时候如果发现$A$可以对角化的时候，可以采用该方法。

如果$k\to\infty$，则$A^k\to 0$（趋于稳定）的条件是什么？从$S\Lambda^kS^{-1}$易得，$|\lambda_i|<1$。再次强调，所有运算的前提是矩阵$A$存在$n$个线性无关的特征向量。如果没有$n$个线性无关的特征向量，则矩阵就不能对角化。

```ad-note
title: 一个结论
同时我们也可以发现$A^n$的性质:如果$A$可以对角化，那么$A^n$的特征值为$\lambda^n$，而特征向量不会变化

```

# Section2.差分方程的应用($A^n$x的计算)

题设:
$u_1=Au_0$，$u_{k+1}=Au_k$，要求解出$u_k$.

我们关注的是可以对角化的矩阵，因此$S$是满秩的，也就是说$n$维向量可以用$S$的列向量表示，即用$A$的特征向量来表示。
$$u_0=c_1x_1+c_2x_2+\cdots+c_nx_n=\Bigg[x_1x_2\cdots x_n\Bigg]\begin{bmatrix}c_1\\c_2\\\vdots\\c_n\end{bmatrix}=Sc$$
我们又知道:$A=S\Lambda S^{-1}$,那么$Au_0=A=S\Lambda S^{-1}Sc=S\lambda c=\lambda Sc$,相应的$A^2u_0=\lambda^2Sc$,则$A^ku_0=\lambda^kSc$

那么我们在考虑$A^{100}u_0$的时候可以考虑两种表达方式，都是一样的：
$$SA^{100}u_0=\Lambda^{100}c=c_1\lambda_1^{100}x_1+c_2\lambda_2^{100}x_2+\cdots+c_n\lambda_n^{100}x_n$$
```ad-tip
这里最重要的一点就是知道$u_0$可以拆分成特征向量的线性组合。

```

# Section3.基于差分方程的实例(斐波那契数列)

$0,1,1,2,3,5,8,13,\cdots,F_{100}=?$，我们要求第一百项的公式，并观察这个数列是如何增长的。可以想象这个数列并不是稳定数列，因此无论如何该矩阵的特征值并不都小于一，这样才能保持增长。而他的增长速度，则有特征值来决定。

我们现在已知的条件仅有:$$F_{k+2}=F_{k+1}+F_{k}$$
但是这个是二阶差分方程，我们需要两个方程来解决，因此也写不出$u_{k+1}=Au_{k}$的形式，因为$A$是方阵，而列向量肯定是二维的(因为方程的自变量肯定有至少两个$F_{k+1},F_{k}$)，我们需要添加条件来构造。

## Subsection1.构造矩阵形式

使用一个**小技巧**，令$u_{k}=\begin{bmatrix}F_{k+1}\\F_{k}\end{bmatrix}$，再追加一个方程组成方程组：$\begin{cases}F_{k+2}&=F_{k+1}+F_{k}\\F_{k+1}&=F_{k+1}\end{cases}$，再把方程组用矩阵表达得到$\begin{bmatrix}F_{k+2}\\F_{k+1}\end{bmatrix}=\begin{bmatrix}1&1\\1&0\end{bmatrix}\begin{bmatrix}F_{k+1}\\F_{k}\end{bmatrix}$，于是我们得到了$u_{k+1}=Au_{k}, A=\begin{bmatrix}1&1\\1&0\end{bmatrix}$。我们把二阶标量方程（second-order scalar problem）转化为一阶向量方程组（first-order system）。

## Subsection2.求解特征值 

我们的矩阵$A=\begin{bmatrix}1&1\\1&0\end{bmatrix}$是一个对称矩阵，==所以它的特征值将会是实数，且他的特征向量将会互相正交==。因为是二阶，我们可以直接利用迹与行列式解方程组$\begin{cases}\lambda_1+\lambda_2&=1\\\lambda_1\cdot\lambda_2&=-1\end{cases}$,$\begin{cases}\lambda_1=\frac{1}{2}\left(1+\sqrt{5}\right)\approx{1.618}\\\lambda_2=\frac{1}{2}\left(1-\sqrt{5}\right)\approx{-0.618}\end{cases}$

我们先来观察这个数列是如何增长的，数列增长由什么来控制？——特征值。哪一个特征值起决定性作用？——较大的一个。

$F_{100}=c_1\left(\frac{1+\sqrt{5}}{2}\right)^{100}+c_2\left(\frac{1-\sqrt{5}}{2}\right)^{100}\approx c_1\left(\frac{1+\sqrt{5}}{2}\right)^{100}$，由于$-0.618$在幂增长中趋近于$0$，所以近似的忽略该项，剩下较大的项，我们可以说数量增长的速度大约是$1.618$。可以看出，这种问题与求解$Ax=b$不同，这是一个动态的问题，$A$的幂在不停的增长，而问题的关键就是这些特征值。

## Subsection3.求解特征向量

求解特征向量，$A-\lambda I=\begin{bmatrix}1-\lambda&1\\1&1-\lambda\end{bmatrix}$，因为有根式且矩阵只有二阶，我们直接观察$\begin{bmatrix}1-\lambda&1\\1&1-\lambda\end{bmatrix}\begin{bmatrix}?\\?\end{bmatrix}=0$，由于$\lambda^2-\lambda-1=0$，则其特征向量为$\begin{bmatrix}\lambda\\1\end{bmatrix}$，即$x_1=\begin{bmatrix}\lambda_1\\1\end{bmatrix}, x_2=\begin{bmatrix}\lambda_2\\1\end{bmatrix}$。

## Subsection4.把$u_0$改写为特征向量的线性组合

最后，计算初始项$u_0=\begin{bmatrix}F_1\\F_0\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix}$，现在将初始项用特征向量表示出来$\begin{bmatrix}1\\0\end{bmatrix}=c_1x_1+c_2x_2$，计算系数得$c_1=\frac{\sqrt{5}}{5}, c_2=-\frac{\sqrt{5}}{5}$。

## Subsection5.计算$u_k=A^ku_0$

因为要求解$F_{100}$，而$u_{k}=\begin{bmatrix}F_{k+1}\\F_{k}\end{bmatrix}$,因此求解$u_{100}$或者$u_{99}$都可以，
我们考虑后者:$u_{99}=S\Lambda^{99}c$，且$c=\begin{bmatrix}c_1\\c_2\end{bmatrix}$有:
$$\begin{aligned}u_{99}=&\begin{bmatrix}F_{100}\\F_{99}\end{bmatrix}=\begin{bmatrix}\frac{1+\sqrt{5}}{2}&\frac{1-\sqrt{5}}{2}\\1&1\end{bmatrix}\begin{bmatrix}\left(\frac{1+\sqrt{5}}{2}\right)^{99}&0\\0&\left(\frac{1-\sqrt{5}}{2}\right)^{99}\end{bmatrix}\begin{bmatrix}\frac{\sqrt{5}}{5}\\-\frac{\sqrt{5}}{5}\end{bmatrix}\\&=\begin{bmatrix}\lambda_1^{99}c_1x_1+\lambda_2^{99}c_2x_2\end{bmatrix}\stackrel{x_1和x_2第一行为两个特征值}{=}\begin{bmatrix}c_1\lambda_1^{100}+c_2\lambda_2^{100}\\c_1\lambda_1^{99}+c_2\lambda_2^{99}\end{bmatrix}\end{aligned}$$
因此:$F_{100}=c_1\lambda_1^{100}+c_2\lambda_2^{100}$，通解为$u_k=c_1\lambda^kx_1+c_2\lambda^kx_2$。

