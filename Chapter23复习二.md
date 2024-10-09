---
author: joe
rating: 8
content: 本文回顾了从复习一到复习二的诸多内容，例如标准正交矩阵，投影矩阵的形式，基于其的施密特正交化，行列式的求法，特征值与特征向量，差分方程。并用题目对其做了一个串联与记忆。
coverage: p25
tags:
  - math
  - linear_algebra
  - class
time: 2023-09-07
---
# Section1.对于复习一到复习二中间的总结

* 我们学习了正交性，有矩阵$Q=\Bigg[q_1\ q_2\ \cdots\ q_n\Bigg]$，若其列向量相互正交，则该矩阵满足$Q^TQ=I$。
* 进一步研究投影，我们了解了Gram-Schmidt正交化法，核心思想是求法向量，即从原向量中减去投影向量$E=b-P, P=Ax=\frac{A^Tb}{A^TA}\cdot A$。
* 接着学习了行列式，根据行列式的前三条性质，我们拓展出了性质4-10。
* 我们继续推导出了一个利用代数余子式求行列式的公式。
* 又利用代数余子式推导出了一个求逆矩阵的公式。
* 接下来我们学习了特征值与特征向量的意义：$Ax=\lambda x$，进而了解了通过$\det(A-\lambda I)=0$求特征值、特征向量的方法。
* 有了特征值与特征向量，我们掌握了通过公式$AS=\Lambda S$对角化矩阵，同时掌握了求矩阵的幂$A^k=S\Lambda^kS^{-1}$。

# Section2.例题

## T1.求$a=\begin{bmatrix}2\\1\\2\end{bmatrix}$的投影矩阵$P$：

   可知$A^TAx=A^Tb\rightarrow \hat{x}=(A^TA)^{-1}A^Tb \rightarrow \hat{b}=A(A^TA)^{-1}A^Tb$，则$P=A(A^TA)^{-1}A^T$,因此:
   $$\begin{aligned}P=P=A&(A^TA)^{-1}A^T\\A^TA=9\ AA^T&=\begin{bmatrix}4&2&4\\2&1&2\\4&2&4\end{bmatrix}\\P=\frac{1}{9}&\begin{bmatrix}4&2&4\\2&1&2\\4&2&4\end{bmatrix}\end{aligned}$$
   ### 小问:
   1. 求$P$矩阵的特征值:
      观察矩阵易知矩阵奇异，且为秩一矩阵，则其零空间为$2$维，所以由$Px=0x$得出矩阵的两个特征值为$\lambda_1=\lambda_2=0$；而从矩阵的迹得知$trace(P)=1=\lambda_1+\lambda_2+\lambda_3=0+0+1$，则第三个特征值为$\lambda_3=1$。

```ad-note
title: 注意对秩为1的矩阵的特征值求法
秩为1的矩阵特征值求法，显然$0$是其中一个，并且我们可以观察$P$的零空间维数，从而判断$0$的可能维数，接着用矩阵的迹来判断另外的特征值

```

2. 求$\lambda_3=1$的特征向量：
	  当然可以直接代入计算，但是我们==对于特殊矩阵可以采用特殊方法==，由$Px=x$我们知道经其意义为，$x$过矩阵$P$变换后不变，又有$P$是向量$a$的投影矩阵，所以任何向量经过$P$变换都会落在$a$的列空间中，则只有已经在$a$的列空间中的向量经过$P$的变换后保持不变，即其特征向量为$x=a=\begin{bmatrix}2\\1\\2\end{bmatrix}$，也就是$Pa=a$。

3. 有差分方程$u_{k+1}=Pu_k,\ u_0=\begin{bmatrix}9\\9\\0\end{bmatrix}$，求解$u_k$：
	   我们当然可以采用差分方程的基础方法，即$A^ku_0=c_1\lambda_1^kx_1+\cdots+c_n\lambda_n^kx_n$，由于投影矩阵的特征值形式都相当好，因此我们也可以这样做，但是还是那句话，==对于特殊矩阵可以采用特殊方法==，我们做了一次$u_1=Pu_0$，就已经把$u_0$投影在了$a$的列空间中，后面再投影的时候已经不会改变了，因此:$u_k=P^ku_0=Pu_0=\begin{bmatrix}6\\3\\6\end{bmatrix}$。

## T2.将点$(1,4),\ (2,5),\ (3,8)$拟合到一条过零点的直线上：

这实际上就是[[Chapter14基于投影矩阵的最小二乘法原理，标准正交向量组#Section2.一个实例(最小二乘法)|最小二乘法]]，设直线为$y=Dt$，写成矩阵形式为$\begin{bmatrix}1\\2\\3\end{bmatrix}D=\begin{bmatrix}4\\5\\8\end{bmatrix}$，即$AD=b$，很明显$D$不存在。
因此我们为了解出该方程，即求最优解，有:
$$A^TA\hat D=A^Tb\rightarrow 14D=38,\ \hat D=\frac{38}{14}\rightarrow y=\frac{38}{14}t$$

## T3.求$a_1=\begin{bmatrix}1\\2\\3\end{bmatrix}\ a_2=\begin{bmatrix}1\\1\\1\end{bmatrix}$的正交向量:

这其实是对[[Chapter15正交矩阵与Gram-Schmidt正交化#Section2.Gram-Schmidt正交化|施密特正交化]]的考察，因此我们利用施密特正交化:

$$\begin{aligned}\beta_1&=\alpha_1=\begin{bmatrix}1\\2\\3\end{bmatrix}\\\beta_2&=\alpha_1-\frac{\alpha_1^T\alpha_2}{\alpha_1^T\alpha_1}\alpha_1=\begin{bmatrix}\frac{4}{7}\\\frac{1}{7}\\\frac{-2}{7}\end{bmatrix}\\&再做正交化即可\end{aligned}$$
## T4.*有$4\times 4$矩阵$A$，其特征值为$\lambda_1,\lambda_2,\lambda_3,\lambda_4$，则矩阵可逆的条件是什么?

矩阵可逆，则行列式显然不能等于0，因此任意的特征值不能为0。
或者理解为零空间中只有零向量，即$Ax=0x$没有非零解，则零不是矩阵的特征值。

### 小问:

1.  *$\det A^{-1}$是什么*：$\det A^{-1}=\frac{1}{\det A}$，而$\det A=\lambda_1\lambda_2\lambda_3\lambda_4$，所以有$\det A^{-1}=\frac{1}{\lambda_1\lambda_2\lambda_3\lambda_4}$。
2. *$trace(A+I)$的迹是什么*：我们知道$trace(A)=a_{11}+a_{22}+a_{33}+a_{44}=\lambda_1+\lambda_2+\lambda_3+\lambda_4$，所以有$trace(A+I)=a_{11}+1+a_{22}+1+a_{33}+1+a_{44}+1=4+\lambda_1+\lambda_2+\lambda_3+\lambda_4$。

```ad-note
title: 关于逆矩阵的特征值
$Ax=\lambda x\rightarrow A^{-1}Ax=\lambda A^{-1}x\rightarrow A^{-1}x=\frac{1}{\lambda}x$
则$A^{-1}$的特征值为$\frac{1}{\lambda}$。

```

## T5.有矩阵$A_4=\begin{bmatrix}1&1&0&0\\1&1&1&0\\0&1&1&1\\0&0&1&1\end{bmatrix}$，求$D_n=?D_{n-1}+?D_{n-2}$：求递归式的系数:

使用代数余子式将矩阵安第一行展开得$\det A_4=1\cdot\begin{vmatrix}1&1&0\\1&1&1\\0&1&1\end{vmatrix}-1\cdot\begin{vmatrix}1&1&0\\0&1&1\\0&1&1\end{vmatrix}=1\cdot\begin{vmatrix}1&1&0\\1&1&1\\0&1&1\end{vmatrix}-1\cdot\begin{vmatrix}1&1\\1&1\end{vmatrix}=\det A_3-\det A_2$。则可以看出有规律$D_n=D_{n-1}-D_{n-2}, D_1=1, D_2=0$。

### 小问:

使用我们在差分方程中的知识构建方程组$\begin{cases}D_n&=D_{n-1}-D_{n-2}\\D_{n-1}&=D_{n-1}\end{cases}$，用矩阵表达有$\begin{bmatrix}D_n\\D_{n-1}\end{bmatrix}=\begin{bmatrix}1&-1\\1&0\end{bmatrix}\begin{bmatrix}D_{n-1}\\D_{n-2}\end{bmatrix}$。计算系数矩阵$A_c$的特征值，$\begin{vmatrix}1-\lambda&1\\1&-\lambda\end{vmatrix}=\lambda^2-\lambda+1=0$，解得$\lambda_1=\frac{1+\sqrt{3}i}{2},\lambda_2=\frac{1-\sqrt{3}i}{2}$，特征值为一对共轭复数。

要判断递归式是否收敛，需要计算特征值的模，即实部平方与虚部平方之和$\frac{1}{4}+\frac{3}{4}=1$。它们是位于单位圆$e^{i\theta}$上的点，即$\cos\theta+i\sin\theta$，从本例中可以计算出$\theta=60^\circ$，也就是可以将特征值写作$\lambda_1=e^{i\pi/3},\lambda_2=e^{-i\pi/3}$。注意，从复平面单位圆上可以看出，这些特征值的六次方将等于1：$e^{2\pi i}=e^{2\pi i}=1$。继续深入观察这一特性对矩阵的影响，$\lambda_1^6=\lambda^6=1$，则对系数矩阵有$A_c^6=I$。则系数矩阵$A_c$服从周期变化，既不发散也不收敛。 

## T6.这样一类矩阵$A_4=\begin{bmatrix}0&1&0&0\\1&0&2&0\\0&2&0&3\\0&0&3&0\end{bmatrix}$，求投影到$A_3$列空间的投影矩阵。

有$A_3=\begin{bmatrix}0&1&0\\1&0&2\\0&2&0\end{bmatrix}$

通常我们会利用$P=A\left(A^TA\right)A^T$来求，但是我们注意到很重要的一点，==如果A是列满秩的话，那么列空间就是满维的空间，投影就会保持不变==，所以按行展开求行列式$\det A_4=-1\cdot-1\cdot-3\cdot-3=9$，所以矩阵可逆，则$P=I$。

求$A_3$的特征值及特征向量：
$\left|A_3-\lambda I\right|=\begin{vmatrix}-\lambda&1&0\\1&-\lambda&2\\0&2&-\lambda\end{vmatrix}=-\lambda^3+5\lambda=0$，解得$\lambda_1=0,\lambda_2=\sqrt 5,\lambda_3=-\sqrt 5$。

我们可以猜测这一类矩阵的规律：奇数阶奇异，偶数阶可逆。

