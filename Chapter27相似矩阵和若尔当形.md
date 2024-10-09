---
author: joe
rating: 8
content: 本章介绍了相似的基本概念，重点介绍并说明了什么情形下两个矩阵能够相似，分为了能对角化前提下与不能对角化前提下，后者引入了Jordan标准型，同时还介绍了相似的性质。
coverage: p29
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-08
---
# Section1.相似矩阵

## Subsection1.定义

矩阵$A,\ B$对于某矩阵$M$满足$B=M^{-1}AM$时，称$A,\ B$互为相似矩阵。

## Subsection2.相似矩阵的存在情形

### 情形1.能够相似对角化时

矩阵何时才能相似对角化？
![[Chapter20矩阵对角化，差分方程(A的乘方的简便计算)#^fb51de|情况1]]![[Chapter20矩阵对角化，差分方程(A的乘方的简便计算)#^45af7c|情况2]]
此时$S^{-1}AS=\Lambda$，则有$A$相似于$\Lambda$。

==如果两个可以相似对角化的矩阵的特征值一样，那他们就相似==，
证明:
$$\begin{aligned}A=S\lambda S^{-1}\ \ B&=C\lambda C^{-1}\\\lambda=S^{-1}AS\ \ \therefore B=&CS^{-1}ASC^{-1}\\令M=(C&S)^{-1}即可\end{aligned}$$
### 情形2.不能相似对角化的时候

例如$\begin{bmatrix}4&1\\0&4\end{bmatrix}$，$\begin{bmatrix}4&1\\0&4\end{bmatrix},\ \begin{bmatrix}5&1\\-1&3\end{bmatrix},\ \begin{bmatrix}4&0\\17&4\end{bmatrix}$，这些矩阵均满足$trace(A)=8,\ \det A=16$。

```ad-note
title: 相似矩阵的共性
相似矩阵有相同的特征值！

证明:
有$Ax=\lambda x,\ B=M^{-1}AM$，第一个式子化为$AMM^{-1}x=\lambda x$，接着两边同时左乘$M^{-1}$得$M^{-1}AMM^{-1}x=\lambda M^{-1}x$，进行适当的分组得$\left(M^{-1}AM\right)M^{-1}x=\lambda M^{-1}x$即$BM^{-1}x=\lambda M^{-1}x$,因此我们发现$B$的特征值没变，特征向量为$M^{-1}x$.

```

我们关注一个特别糟糕的矩阵:
$$\begin{bmatrix}0&1&0&0\\0&0&1&0\\0&0&0&0\\0&0&0&0\end{bmatrix}$$其特征值为四个零。很明显矩阵的秩为$2$，所以其零空间的维数为$4-2=2$，即该矩阵有两个特征向量。可以发现该矩阵在主对角线的上方有两个$1$，在对角线上每增加一个$1$，特征向量个个数就减少一个。

再提出一个矩阵:
$$\begin{bmatrix}0&1&0&0\\0&0&0&0\\0&0&0&1\\0&0&0&0\end{bmatrix}$$
他的特征值也是四个0，特征向量也是2个，从特征向量的数目看来这两个矩阵是相似的，其实不然。

# Section2.$Jordan$形

基于上述两个矩阵，若尔当认为第一个矩阵是由一个$3\times 3$的块与一个$1\times 1$的块组成的 $\left[\begin{array}{ccc|c}0&1&0&0\\0&0&1&0\\0&0&0&0\\\hline0&0&0&0\end{array}\right]$，而第二个矩阵是由两个$2\times 2$矩阵组成的$\left[\begin{array}{cc|cc}0&1&0&0\\0&0&0&0\\\hline0&0&0&1\\0&0&0&0\end{array}\right]$，这些分块(指的是左上角和右下角的那两个)被称为若尔当块。

那么由于第一个矩阵的若尔当块为$J_1和J_3$，而第二个矩阵的若尔当块为两个$J_2$，因此二者不相似。

## Subsection1.概念

若尔当块的定义型为$J_i=\begin{bmatrix}\lambda_i&1&&\cdots&\\&\lambda_i&1&\cdots&\\&&\lambda_i&\cdots&\\\vdots&\vdots&\vdots&\ddots&\\&&&&\lambda_i\end{bmatrix}$，它的对角线上只为同一个数因此特征值只有1个，仅有一个特征向量。

若尔当矩阵的定义为$J=\left[\begin{array}{c|c|c|c}J_1&&&\\\hline&J_2&&\\\hline&&\ddots&\\\hline&&&J_d\end{array}\right]$，若尔当块的个数即为矩阵特征值的个数。

## Subsection2.结论

Jordan证明了，任何方阵都能化为$Jordan$标准型，$M$的若尔当标准型可以写成$P^{-1}MP= J$，因此如果两个矩阵的$Jordan$标准型不同，则说明二者不相似。

具体的$Jordan$化过程比较复杂，未作具体的阐述。

但是，==无论两个矩阵是否能够进行相似对角化，只要特征值不一样，那二者就不相似==。

另外，==能够相似对角化的矩阵只能和能相似对角化的矩阵相似，不能的和不能的相似==

对于第二个做一个简单的证明:
$$A=S\lambda S^{-1}\ \ \ B=M^{-1}AM=M^{-1}S\lambda S^{-1}M=(S^{-1}M)^{-1}\lambda S^{-1}M$$
因此$B$可以相似对角化。

当然在矩阵情况良好的前提下，$Jordan$标准型也适用，$n$阶矩阵将有$n$个不同的特征值，它的若尔当矩阵就是$\Lambda$，共$n$个特征向量，有$n$个若尔当块。