---
author: joe
rating: 8
content: 关于AX=b的解的情况讨论，矩阵的秩与列满秩行满秩不同关系下方程有什么不同的解情况
coverage: p8
tags:
  - linear_algebra
  - math
  - class
---

# Section1.AX=b什么时候有解？

我们针对同样的系数矩阵[[Chapter05AX=0的具体算法，秩与特解的得到#Subsection1.关于$Ax=0$与$Ux=0$|A]]，改变右侧的目标向量，来考虑是否有解.

$$A=\begin{bmatrix}1&2&2&2\\2&4&6&8\\3&6&8&10\end{bmatrix}$$
而我们考虑:
$$b=\begin{bmatrix}b_1\\b_2\\b_3\end{bmatrix}$$
我们仍然利用消元法，同时我们考虑增广矩阵,有:
$$[U\ |\ b\ ]=\begin{bmatrix}1&2&2&2|&b_1\\0&0&2&4|&b_2-2b_1\\0&0&0&0|&b_3-b_1-b_2\end{bmatrix}$$
那么这里我们观察到，显然最后一行的右侧向量$b_3-b_1-b_2=0$，才能有解.

我们其实在[[Chapter04向量空间#Subsection1.列空间|列空间]]中已经提到过，这里如果$AX=b$要有解，需要$b$在$A$的列空间中，这是一种说法；
根据上面这个例子，我们有另外一种说法，就是:系数矩阵做了某些行变换得到零行，那么$b$的分量做同样的变化也要得到0.

# Section2.Ax=b解的结构:特解＋通解

我们这里让$b=\begin{bmatrix}1\\5\\6\end{bmatrix}$,满足上面的可解性条件$b_3-b_2-b_1=0$，我们会得到两个方程:

$$\begin{aligned}&x_1+2x_2+2x_3+2x_4=1\\&2x_3+4x_4=3\end{aligned}$$

这个时候我们令自由变量$\begin{bmatrix}x_2\\x_4\end{bmatrix}=\begin{bmatrix}0\\0\end{bmatrix}$,那么我们可以得到一个方程的解为$$X_{particular}=\begin{bmatrix}x_1\\x_2\\x_3\\x_4\end{bmatrix}=\begin{bmatrix}-1\\0\\ \frac{3}{2}\\0\end{bmatrix}$$,这个我们称之为方程的特解. 

我们在[[Chapter05AX=0的具体算法，秩与特解的得到]]中已经知道了$AX=0$的解，即

$$X_{null}=k_1\begin{bmatrix}-2\\0\\1\\0\end{bmatrix}+k_2\begin{bmatrix}2\\-1\\0\\1\end{bmatrix}$$
那么我们这里就认为$AX=b$的解就是

$$X=X_{null}+X_{paticular}$$ ^e0e000
## 向量空间的形式:

注意到$X_{null}=k_1\begin{bmatrix}-2\\0\\1\\0\end{bmatrix}+k_2\begin{bmatrix}2\\-1\\0\\1\end{bmatrix}$，这是一个位于$R^4$的$R^2$平面(因为只有$k_1和k_2$两个变量，那么只能构成一个平面)，而$X_{particular}=\begin{bmatrix}x_1\\x_2\\x_3\\x_4\end{bmatrix}=\begin{bmatrix}-1\\0\\ \frac{3}{2}\\0\end{bmatrix}$是一个直线向量，所以相当于一个平面上的任意向量加一个直线向量.

## prove:为什么AX=b的结构是特解＋零空间的通解?

一个简单的证明:

$$\begin{aligned}\because AX_{paticular}=b\ \ AX_{null}=0\\
\therefore A(X_{paticular}+X_{null})=b\end{aligned}$$
因此我们可以知道$AX=b$的解为
$$X=X_{p}+X_{n}$$
## addition:为什么我们取$\begin{bmatrix}x_2\\x_4\end{bmatrix}=\begin{bmatrix}0\\0\end{bmatrix}$来求特解?

因为我们取其他的值，会发现这个与$X_{null}$重复(这个证明很容易，因为我们得到$X_{null}$是通过利用$\begin{bmatrix}x_2\\x_4\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix}$来得到的，我们取其他值只不过是$k$倍的罢了.)

# Section2.行满秩，列满秩，矩阵秩，解的形式之间的关系.

## 列满秩，$n=r$

^017718

$$列秩n=r\rightarrow变量数目＝有用的方程数\rightarrow没有自由变量\rightarrow零空间只有零向量\rightarrow 如果b符合要求，AX=b只有一个解，即特解$$
## 行满秩，$m=r$
$$行秩m=r\rightarrow 每一行一定会有主元\rightarrow对于b没有要求\rightarrow AX=b一定会有解$$
## 列满秩且行满秩
$$m=n=r意味着一定是方阵\rightarrow AX=b一定有解，为特解$$


| r=m=n                             | r=m<n                               | r=n<m                                   | r<n,r<m                                  |
| --------------------------------- | ----------------------------------- | --------------------------------------- | ---------------------------------------- |
| $R=\begin{bmatrix}I\end{bmatrix}$ | $R=\begin{bmatrix}I&F\end{bmatrix}$ | $R=\begin{bmatrix}I\\0\end{bmatrix}$    | $R=\begin{bmatrix}I&F\\0&0\end{bmatrix}$ |
|  存在唯一解，即特解                | 一定有解，且有无穷多解              | 对零行的b有要求，如果符合要求则为唯一解 | 0解(b不符合要求)或无穷解                 |

^618d69
