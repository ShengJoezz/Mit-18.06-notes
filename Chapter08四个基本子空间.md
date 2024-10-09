---
author: joe
rating: 6
content: 这一章具体介绍了系数矩阵A的四个基本子空间
coverage: p10
tags:
  - linear_algebra
  - math
  - class
---
# Section1.引入

在[[Chapter04向量空间]]中我们介绍了列空间与零空间两个基本的子空间，我们考虑这两个空间的因素也是基于$AX=b$的求解需求来考虑的，而右乘是列变换，因此我们只考虑了列空间与零空间，事实上，我们不仅仅是利用右乘来考虑问题，行变换也是我们要利用的基本变换，因此在这里我们引入行空间与左零空间的说法。

即四个基本子空间为:
* 列空间($C(A)$)
* 零空间($N(A)$)
* 行空间($C(A^T)$)
* 左零空间($N(A^T)$)

## 左零空间是什么？

我们考虑行变换的零空间$N(A^T)$，需要考虑:$A^Ty=0$，但是我们往往是考察$A$而非$A^T$的性质，因此我们同时取转置，可以得到:

$$y^TA=0$$

这个$y^T$在左侧，因此我们称之为左零空间。

==by the way,我觉得把行空间考虑成为$A^T$的列空间是一件很coooooool的事情

# Section2.四个基本子空间的维数

## Subsection1.列空间和行空间的维数

我们知道[[Chapter05AX=0的具体算法，秩与特解的得到#^7af8fa|列秩=行秩=矩阵秩=r]]:

$$A的主元个数=A的列向量的线性无关的数目(列的维数)=A的行向量的线性无关的数目(行的维数)=A^T的主元个数=r$$

那么我们可以知道$dimension(行空间)=dimension(列空间)=r$

## Subsection2.零空间和左零空间

关于$dimension(N(A))$,我们在[[Chapter05AX=0的具体算法，秩与特解的得到#Subsection1.方程解向量的数量]]中提到，解向量的数目正是自由变量的数目，因此:

$$dimension(N(A))=n-r$$
而$dimension(N(A^T))$是类似的，
$$dimension(N(A^T))=m-r$$

## Subsection3.一个发现

$$dimension(C(A))+dimension(N(A))=n$$
$$dimension(C(A^T))+dimension(N(A^T))=m$$

事实上这也是很好理解的，$A$的列空间的基向量数目依赖于主元的个数，而零空间的基向量数目依赖于自由变量的数目，相加刚好就是列数.

# Section3.四个基本子空间的基的求法

## Subsection1.列空间$C(A)$

利用消元法把$A$变成$RREF$之后，主元所在的列即为$C(A)$.

## Subsection2.零空间$N(A)$

根据[[Chapter05AX=0的具体算法，秩与特解的得到#^bbe1c6|Ax=0的求解]]我们可以顺利得到$AX=0$的解向量.

## Subsection3.行空间$C(A^T)$

==很有意思的想法.

在[[#Subsection1.列空间$C(A)$|列空间的基求法]]中我们提到利用行变换把$A$变成$RREF$，主元所在的列即为$C(A)$.

但是我们注意到一点，我们在做行变换的时候，行空间始终没有变化，得到某些零行实际就是把线性相关的行消除了，得到的非0行实际上就是线性独立的.

==amazing！我们利用行变换把$A$变成$RREF$的过程竟然能够把$C(A)$和$C(A^T)$全部得到。

## Subsection4.左零空间$C(A^T)$

^129864

我们注意到$RREF$的样式:

$$RREF=\begin{bmatrix}I&F\\0&0\end{bmatrix}$$
再关注所需要求$C(A^T)$的格式:

$$y^TA=0$$
我们注意到我们在得到$RREF$的过程中，我们利用了行变换进行消元，根据[[Chapter02矩阵的求逆与A=LU分解法#^5c9a4c|Gauss-Jordan方法]]，提示我们利用初等矩阵来实现行变换，即有:
$$E\times A=RREF$$
我们发现$RREF$的最下几行是$0$，那么说明$E$中存储了某个组合，使得$A$中的行向量组合为$0$

那么我们考虑利用[[Chapter02矩阵的求逆与A=LU分解法#^5c9a4c|Gauss-Jordan方法]]，利用$\begin{bmatrix}A&I\end{bmatrix}$,来进行同步变化，得到$\begin{bmatrix}RREF&E\end{bmatrix}$，观察相关的行为0对应的$E$,这个就是$A$的左零空间的基。

