---
author: joe
rating: 8
content: 本章对前面的内容利用题目做了集中的回顾，尤其是对于方程与系数矩阵秩的关系，这里我做了更细致的认识。
coverage: p35
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-10
---
# Section1.试题们

## T1.方程与系数矩阵秩的关系

已知$m\times n$矩阵$A$，有$Ax=\begin{bmatrix}1\\0\\0\end{bmatrix}$无解；$Ax=\begin{bmatrix}0\\1\\0\end{bmatrix}$仅有唯一解，求关于$m,n,rank(A)$的信息。

```ad-note
title: 矩阵与秩的关系
$n$是用来控制零空间/$Ax=0$的解情况：
1.如果$r>n$，则$Ax=0$有无穷解，$Ax=b$未必有解，如果有解则有无穷解。
2.如果$r=n$，则$Ax=0$只有$x=0$这一解,$Ax=b$未必有解，如果有解则有唯一解。

$m$是用来控制$Ax=b$的解情况:
1.如果$r>m$，则$Ax=b$未必有解，因为$A$消元会产生零行，相应的$b$也需要是0
2.如果$r=m$，则$Ax=b$一定有解，因为$A$的主元每一行都有，则一定有解。

```

1. 由于$Ax=\begin{bmatrix}1\\0\\0\end{bmatrix}$无解，那么$m$是不能等于$r$的，则$m<r$；
又因为 $Ax=\begin{bmatrix}0\\1\\0\end{bmatrix}$仅有唯一解，那么$Ax=0$只有零解，因此$n=r=3$，故$m<n=r=3$。

写出一个矩阵$A$的特例：$A=\begin{bmatrix}0&0\\1&0\\0&1\end{bmatrix}$。

2. $\det A^TA\stackrel{?}{=}\det AA^T$
   答案是否定的，因为$A^TA的秩情况与A的列情况一样$，而$AA^T的秩情况与A的行情况一样$，因此$r(A^TA)=3,r(AA^T)=2$，则第一个行列式不等于0，而第二个等于0.（但是对于方阵，$\det AB=\det BA$恒成立。）
   
3. $A^TA$可逆吗？
   由2可知，$r(A^TA)=3$，因此可逆。

4. $AA^T$正定吗？
   显然是不正定的，因为$r(AA^T)=2$,所以会有$\lambda=0$，故不正定，而是半正定矩阵。

5. 求证$A^Ty=c$至少有一个解
   可知$A^T=n\times m$的矩阵，因此行满秩，一定有解，且$n-r>0$，因此零空间也不为0，因此该方程至少有一个解。

## T2.一道简单的题目
设$A=\Bigg[v_1\ v_2\ v_3\Bigg]$，对于$Ax=v_1-v_2+v_3$，求$x$。
 按列计算矩阵相乘，有$x=\begin{bmatrix}1\\-1\\1\end{bmatrix}$。
 
 若$Ax=v_1-v_2+v_3=0$，则解是唯一的吗？为什么。
 如果解释唯一的，则零空间中只有零向量，而在此例中$x=\begin{bmatrix}1\\-1\\1\end{bmatrix}$就在零空间中，所以解不唯一，$kx$也是解。
 
  若$v_1,v_2,v_3$是标准正交向量，那么怎样的线性组合$c_1v_1+c_2v_2$能够最接近$v_3$？此问是考察投影概念，由于是正交向量，所以只有$0$向量最接近$v_3$。

## T3.马尔可夫矩阵

矩阵$A=\begin{bmatrix}.2&.4&.3\\.4&.2&.3\\.4&.4&.4\end{bmatrix}$，求稳态。


这是个马尔科夫矩阵，前两之和为第三列的两倍，奇异矩阵总有一个特征值为$0$，而马尔科夫矩阵总有一个特征值为$1$，剩下一个特征值从矩阵的迹得知为$-.2$。
再看马尔科夫过程，设从$u(0)$开始，$u_k=A^ku_0, u_0=\begin{bmatrix}0\\10\\0\end{bmatrix}$。先代入特征值$\lambda_1=0,\ \lambda_2=1,\ \lambda_3=-.2$查看稳态$u_k=c_1\lambda_1^kx_1+c_2\lambda_2^kx_2+c_3\lambda_3^kx_3$，当$k\to\infty$，第一项与第三项都会消失，剩下$u_\infty=c_2x_2$。
  
到这里我们只需求出$\lambda_2$对应的特征向量即可，带入特征值求解$(A-I)x=0$，有$\begin{bmatrix}-.8&.4&.3\\.4&-.8&.3\\.4&.4&-.6\end{bmatrix}\begin{bmatrix}?\\?\\?\end{bmatrix}=\begin{bmatrix}0\\0\\0\end{bmatrix}$，可以消元得，也可以直接观察得到$x_2=\begin{bmatrix}3\\3\\4\end{bmatrix}$。 

剩下就是求$c_2$了，可以通过$u_0$一一解出每个系数，但是这就需要解出每一个特征值。另一种方法，我们可以通过马尔科夫矩阵的特性知道，对于马尔科夫过程的每一个$u_k$都有其分量之和与初始值分量之和相等，所以对于$x_2=\begin{bmatrix}3\\3\\4\end{bmatrix}$，有$c_2=1$。所以最终结果是$u_\infty=\begin{bmatrix}3\\3\\4\end{bmatrix}$。

## T4.二阶方阵(已知特征值和特征向量怎么得到原矩阵)

 1. 求投影在直线$a=\begin{bmatrix}4\\-3\end{bmatrix}$上的投影矩阵：应为$P=\frac{aa^T}{a^Ta}$。

2. 已知特征值$\lambda_1=2,\ x_1=\begin{bmatrix}1\\2\end{bmatrix}\quad \lambda_2=3,\ x_2=\begin{bmatrix}2\\1\end{bmatrix}$求原矩阵$A$：从对角化公式得$A=S\Lambda S^{-1}=\begin{bmatrix}1&2\\2&1\end{bmatrix}\begin{bmatrix}0&0\\0&3\end{bmatrix}\begin{bmatrix}1&2\\2&1\end{bmatrix}^{-1}$，解之即可。
3. $A$是一个实矩阵，且对任意矩阵$B$，$A$都不能分解成$A=B^TB$，给出$A$的一个例子：我们知道$B^TB$是对称的，所以给出一个非对称矩阵即可。
4. 非对称矩阵的特征向量可以是正交的吗？
   可以的，反对称矩阵，因为满足$AA^T=A^TA$而同样具有正交的特征向量，所以有$A=\begin{bmatrix}0&1\\-1&0\end{bmatrix}$或旋转矩阵$\begin{bmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{bmatrix}$，这些矩阵都具有复数域上的正交特征向量组。

## T5.最小二乘法问题

最小二乘问题，因为时间的关系直接写出计算式和答案，$\begin{bmatrix}1&0\\1&1\\1&2\end{bmatrix}\begin{bmatrix}C\\D\end{bmatrix}=\begin{bmatrix}3\\4\\1\end{bmatrix}(Ax=b)$，解得$\begin{bmatrix}\hat C\\\hat D\end{bmatrix}=\begin{bmatrix}\frac{11}{3}\\-1\end{bmatrix}$。

求投影后的向量$p$：向量$p$就是向量$b$在矩阵$A$列空间中的投影，所以$p=\begin{bmatrix}p_1\\p_2\\p_3\end{bmatrix}=\begin{bmatrix}1&0\\1&1\\1&2\end{bmatrix}\begin{bmatrix}\hat C\\\hat D\end{bmatrix}$。

求拟合直线的图像：$x=0,1,2$时$y=p_1,p_2,p_2$所在的直线的图像，$y=\hat C+\hat Dx$即$y=\frac{11}{3}-x$。

![[Pasted image 20230910171827.png|379]]

求一个向量$b$使得最小二乘求得的$\begin{bmatrix}\hat C\\\hat D\end{bmatrix}=\begin{bmatrix}0\\0\end{bmatrix}$：
我们知道最小二乘求出的向量$\begin{bmatrix}\hat C\\\hat D\end{bmatrix}$使得$A$列向量的线性组合最接近$b$向量（即$b$在$A$列空间中的投影），如果这个线性组合为$0$向量（即投影为$0$），则$b$向量与$A$的列空间正交，所以可以取$b=\begin{bmatrix}1\\-2\\1\end{bmatrix}$同时正交于$A$的两个列向量。

