---
author: joe
rating: 8
content: 本章介绍了实对称的三个性质并提到了一定的证明，同时介绍了其性质更好的一个子类正定矩阵的四个性质。
coverage: p26
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-07
---
# Section1.实对称矩阵

## Subsection1,实对称矩阵的特性

1. 特征值为实数；（对比第二十一讲介绍的旋转矩阵，其特征值为纯虚数。）
2. 特征向量相互正交。（当特征值重复时，特征向量也可以从子空间中选出相互正交正交的向量。）
3. ==主元符号的正负数量与特征值的正负数量相同==。

### 证明:
$$\begin{aligned}Ax=\lambda x\stackrel{取共轭}{\rightarrow}\bar A\bar x=\bar\lambda \bar x\stackrel{A是实数矩阵}{\rightarrow}A\bar x=\bar\lambda \bar x\stackrel{取转置}{\rightarrow}{\bar{x}^TA=\bar{x}^T\bar\lambda}\end{aligned}$$
那么我们得到两个等式:$Ax=\lambda x$与${\bar{x}^TA=\bar{x}^T\bar\lambda}$：
那么$$\begin{aligned}Ax=\lambda x\stackrel{左乘\bar{x}^T}{\rightarrow}&=\bar{x}^TAx=\bar{x}^T\lambda x\\\bar{x}^TA=\bar{x}^T\bar\lambda\stackrel{右乘x}{\rightarrow}&\bar{x}^TAx=\bar{x}^T\bar\lambda x\\\therefore\bar{x}^T\lambda x=&\bar{x}^T\bar\lambda x\\\therefore\lambda=\bar{\lambda}(\bar{x}^Tx\neq 0)&则\lambda虚部为0\end{aligned}$$
这里有一个假设，就是$\bar{x}^Tx\neq 0$，如果$x$是实数向量显然是易证的，但是此时$x$是复数向量，因此我们还要一个证明:

$\bar{x}^Tx=\begin{bmatrix}\bar x_1&\bar x_2&\cdots&\bar x_n\end{bmatrix}\begin{bmatrix}x_1\\x_2\\\vdots\\x_n\end{bmatrix}=\bar x_1x_1+\bar x_2x_2+\cdots+\bar x_nx_n$，设$x_1=a+ib, \bar x_1=a-ib$则$\bar x_1x_1=a^2+b^2$，所以有$\bar{x}^Tx>0$。而$\bar{x}^Tx$就是$x$长度的平方。

```ad-note
title: 拓展(如果A是复数矩阵那何时才能成立呢？)
则需要在第一步取共轭的时候保留$\bar{A}$的共轭形式，因此在进行后面两个相同的量的比较，也需要$A=\bar{A}^T$

```

## Subsection2.实对称矩阵的写法

在通常（可对角化）情况下，一个矩阵可以化为：$A=S\varLambda S^{-1}$,在矩阵对称的情况下，通过性质2可知，由特征向量组成的矩阵$S$中的列向量是相互正交的，此时如果我们把特征向量的长度统一化为$1$，就可以得到一组标准正交的特征向量。则对于对称矩阵有$A=Q\varLambda Q^{-1}$，而对于标准正交矩阵，有$Q=Q^T$，所以对称矩阵可以写为$$A=Q\varLambda Q^T$$
也可以写为:
$$A=Q\varLambda Q^T=\Bigg[q_1\ q_2\ \cdots\ q_n\Bigg]\begin{bmatrix}\lambda_1& &\cdots& \\&\lambda_2&\cdots&\\\vdots&\vdots&\ddots&\vdots\\& &\cdots&\lambda_n\end{bmatrix}\begin{bmatrix}\quad q_1^T\quad\\\quad q_1^T\quad\\\quad \vdots \quad\\\quad q_1^T\quad\end{bmatrix}=\lambda_1q_1q_1^T+\lambda_2q_2q_2^T+\cdots+\lambda_nq_nq_n^T$$

注意这个展开式中的$qq^T$，$q$是单位列向量所以$q^Tq=1$，结合我们在第十五讲所学的投影矩阵的知识有$\frac{qq^T}{q^Tq}=qq^T$是一个投影矩阵，很容易验证其性质，比如平方它会得到$qq^Tqq^T=qq^T$于是多次投影不变等。

**每一个对称矩阵都可以分解为一系列相互正交的投影矩阵。**

# Section2.正定矩阵

## Subsection1.概念

正定矩阵是对称矩阵中性质更好的一个子类(对称矩阵本身的性质就已经很好了)，正定矩阵指特征值均为正数的矩阵（根据上面的性质3有矩阵的主元均为正）。 ^35c749

举个例子，$\begin{bmatrix}5&2\\2&3\end{bmatrix}$，由行列式消元知其主元为$5,\frac{11}{5}$，按一般的方法求特征值有$\begin{vmatrix}5-\lambda&2\\2&3-lambda\end{vmatrix}=\lambda^2-8\lambda+11=0, \lambda=4\pm\sqrt 5$。

## Subsection2.性质

正定矩阵的另一个性质是，所有子行列式为正。对上面的例子有$\begin{vmatrix}5\end{vmatrix}=5, \begin{vmatrix}5&2\\2&3\end{vmatrix}=11$。(这是显然的，因为所有的主元都大于0，那么所有的子行列式包含的主元都会大于0)
```ad-summary
title: 正定矩阵的特性
1.本身就是对称矩阵
2.特征值均大于0
3.主元均大于0(对称矩阵都有这个性质)
4.所有的子行列式均大于0

```

我们看到正定矩阵将早期学习的的消元主元、中期学习的的行列式、后期学习的特征值结合在了一起。