---
author: joe
rating: 8
content: 这一章我们介绍了行列式计算的三种方式，除了我们在上一章节中提到的消元化为三角矩阵后用对角线元素累乘这种方式，我们可以用逆序数参与的排列组合算法以及代数余子式累加算法。
coverage: p19
tags:
  - linear_algebra
  - math
  - class
---
# Section1.行列式公式推导

## Subsection1.行列式的推导

我们已经知道$\begin{vmatrix}a& b\\c& d\end{vmatrix}=ad-bc$，我们可以用很多方法来推导，在[[Chapter16行列式的基本运算法则#性质 8 当且仅当 ${A}$ 是奇异矩阵(不可逆)时,$ det {A}=0$。|消元法改写为上三角矩阵]]我们可以推导，这里为了引进一般维数的行列式公式，我们考虑另外一种证明方式:
$$\begin{vmatrix}a& b\\c& d\end{vmatrix}\stackrel{3.b}{=}\begin{vmatrix}a& 0\\c& d\end{vmatrix}+\begin{vmatrix}0& b\\c& d\end{vmatrix}\stackrel{3.b}{=}\begin{vmatrix}a& 0\\c& 0\end{vmatrix}+\begin{vmatrix}a& 0\\0& d\end{vmatrix}+\begin{vmatrix}0& b\\c& 0\end{vmatrix}+\begin{vmatrix}0& b\\0& d\end{vmatrix}$$
$$其中 \begin{vmatrix}a& 0\\c& 0\end{vmatrix}\stackrel{6, 10}{=}0,\quad\begin{vmatrix}a& 0\\0& d\end{vmatrix}\stackrel{7}{=}ad,\quad \begin{vmatrix}0& b\\c& 0\end{vmatrix}\stackrel{2,7}{=}-bc,\quad \begin{vmatrix}0& b\\0& d\end{vmatrix}\stackrel{6,10}{=}0$$
$$综上，\begin{vmatrix}a& b\\c& d\end{vmatrix}=ad-bc$$
我们注意到我们在拆分的过程中，我们保留某一行/列全为0的行列式是没有意义的，因此我们在拆分的过程中需要注意到每行每列都需要保留元素。我们拆分的最终结果一定会是任意位置的0与任意位置的元素的组合，事实上我们可以一点一点的写开，这里就不做具体的证明。

我们接着拿三阶矩阵来举例子:

$$\begin{align}
\begin{vmatrix}a_{11}& a_{12}& a_{13}\\a_{21}& a_{22}& a_{23}\\a_{31}& a_{32}& a_{33}\end{vmatrix}&=\begin{vmatrix}a_{11}& 0& 0\\0& a_{22}& 0\\0& 0& a_{33}\end{vmatrix}+\begin{vmatrix}a_{11}& 0& 0\\0& 0& a_{23}\\0& a_{32}& 0\end{vmatrix}+\begin{vmatrix}0& a_{12}& 0\\a_{21}& 0& 0\\0& 0& a_{33}\end{vmatrix}\\&+\begin{vmatrix}0& a_{12}& 0\\0& 0& a_{23}\\a_{31}& 0& 0\end{vmatrix}
+\begin{vmatrix}0& 0& a_{13}\\a_{21}& 0& 0\\0& a_{32}& 0\end{vmatrix}+\begin{vmatrix}0& 0& a_{13}\\0& a_{22}& 0\\a_{31}& 0& 0\end{vmatrix}\\
&=a_{11}a_{22}a_{33}-a_{11}a_{23}a_{32}-a_{12}a_{21}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32}-a_{13}a_{22}a_{31}
\end{align}$$
那我们推广一些，考虑$n$阶矩阵，可以知道$n$ 阶行列式可以分解得到 $n!$ 个非零行列式(占据第1行的元素有$n$种选择,占据第2行的元素有$n-1$种选择,以此类推即$n!$)，因此我们可以得到任意的$n$阶矩阵的行列式，为:
$$\det A=\sum_{n!} \pm a_{1\alpha}a_{2\beta}a_{3\gamma}\cdots a_{n\omega}, \quad(\alpha, \beta, \gamma, \cdots,\omega)=P_n^n$$
其中 $\alpha, \beta, \gamma, \cdots, \omega$ 共 $n$ 个数是列标 $1$ 到 $n$ 的某种排列(这意味着不会重复)。==上式即为一般的行列式公式==。

## Subsection2.逆序数

那么我们需要继续考虑，这个正负号怎么取到呢？我们考虑如下的例子:
$$\begin{vmatrix}0& 0& \color{red}1& \color{blue} 1\\0& \color{red} 1& \color{blue} 1& 0\\ \color{red} 1& \color{blue} 1& 0& 0\\\color{blue} 1& 0& 0&  \color{red} 1\end{vmatrix}$$
我们按上面的方法来取元素来相乘，第一行我们容易发现,$\alpha$只能为3或4。选定$\alpha=3$,则$\beta=2,\gamma=1,\omega=4$;选定$\alpha=4$,则$\beta=3,\gamma=2,\omega=1$。我们只能找到两组非零行列式。
- 对于第一组非零行列式(图中蓝色部分),其对应的排列是 $(4,3,2,1)$,变为 $(1,2,3,4)$ 需要两步操作,由行列式性质2可知符号应取 $+$.
- 对于第二组非零行列式(图中红色部分),其对应的排列是 $(3,2,1,4)$,变为 $(1,2,3,4)$ 需要一步操作,由行列式性质2可知符号应取 $-$.
我们会发现，这个正负号的判别依赖于所留存元素的列变换，而这个列变换我们可以抽象为$\alpha, \beta, \gamma, \cdots, \omega$ 这几个数字换成$(1,2,3\cdots n)$的次数，因此我们考虑到==逆序数的概念==:

逆序数就是从左到右遍历当前排列中的每一个数,统计右侧有几个数比自己小,比如对于排列 $(4,3,2,1)$,$4$后面有3个数比它小,$3$后面有2个数比它小,$2$后面有1个数比它小,取其总和$3+2+1=6$即为逆序数，记为$r(j_1,j_2,\cdots,j_n)$。
==逆序数为偶数的,则称排列为偶排列,否则为奇排列,偶排列时非零行列式取 $+$,奇排列时非零行列式取 $-$==。其中的道理在于奇排列做一次交换后即为偶排列,偶排列做一次交换后即为奇排列,并且初始顺序排列 $(1,2,3,\cdots,n)$ 为偶排列。

因此行列式公式可以改写为:
  $$\det A=\sum_{n!} (-1)^{r(j_1,j_2,\cdots,j_n)} a_{1\alpha}a_{2\beta}a_{3\gamma}\cdots a_{n\omega}, \quad(\alpha, \beta, \gamma, \cdots,\omega)=P_n^n$$

# Section2.代数余子式

我们考虑三阶矩阵的行列式，实际上我们已经推导出来它的值:

$$原式 =a_{11}a_{22}a_{33}-a_{11}a_{23}a_{32}-a_{12}a_{21}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32}-a_{13}a_{22}a_{31}$$
我们如果把第一行的元素做同类项合并，会有:
$$\begin{aligned}a_{11}(a_{22}a_{33}-a_{23}a_{32})+a_{12}(-a_{21}a_{33}+a_{23}a_{31})+a_{13}(a_{21}a_{32}-a_{22}a_{31})\\=\begin{vmatrix}a_{11}& 0& 0\\0& a_{22}& a_{23}\\0& a_{32}& a_{33}\end{vmatrix}+\begin{vmatrix}0& a_{12}& 0\\a_{21}& 0& a_{23}\\a_{31}& 0& a_{33}\end{vmatrix}+\begin{vmatrix}0& 0& a_{13}\\a_{21}& a_{22}& 0\\a_{31}& a_{32}& 0\end{vmatrix}\\=a_{11}(\begin{vmatrix} a_{22}& a_{23}\\a_{32}& a_{33}\end{vmatrix})+a_{12}(-\begin{vmatrix} a_{21}& a_{23}\\a_{31}& a_{33}\end{vmatrix})+a_{13}(\begin{vmatrix} a_{21}& a_{22}\\a_{31}& a_{32}\end{vmatrix})\end{aligned}$$
容易发现,3×3 的行列式由 2×2 行列式组成。事实上,$n$阶行列式同样可化为多个 $n-1$ 阶行列式的组合。下面我们正式介绍 $a_{ij}$的==代数余子式 的概念==:
$$C_{ij}=(-1)^{i+j}\cdot det(去掉 i 行和 j 列的一个 n-1 矩阵)$$
因此我们得到新的求解行列式的方式:
$$假设我们选取第一行（i=1），那么行列式 A 沿第一行展开有:\\\det(A)=a_{11}C_{11}+a_{12}C_{12}+...+a_{1n}C_{1n}$$
```ad-note
title: 代数余子式和余子式的区别:
余子式即去掉 $i$行和$j$列的一个$n-1$行列式；

代数余子式需要在余子式的基础上带上符号，$C_{ij}=(-1)^{i+j}\cdot det(去掉 i 行和 j 列的一个 n-1 矩阵)$。

```

# Section3.求解行列式的三种方式

- 消元,将矩阵化为三角矩阵,主元乘积记为行列式的值(最简单)
- 按照行列式公式将行列式完全展开,找到 $n!$ 种非零行列式,计算这些行列式的值的和(最复杂)
  $$\det A=\sum_{n!} (-1)^{r(j_1,j_2,\cdots,j_n)} a_{1\alpha}a_{2\beta}a_{3\gamma}\cdots a_{n\omega}, \quad(\alpha, \beta, \gamma, \cdots,\omega)=P_n^n$$
- 使用代数余子式对行列式进行降阶,展开得到更简单的行列式,然后再求解(介于二者之间)
  $$假设我们选取第一行（i=1），那么行列式 A 沿第一行展开有:\\\det(A)=a_{11}C_{11}+a_{12}C_{12}+...+a_{1n}C_{1n}$$

# Section4.一个有趣的例子

我们注意到一个特殊的矩阵，叫做三对角矩阵，这里我们考虑的是由1组成的:

$$A_{1}=1, A_{2}=\begin{bmatrix}{1} & {1} \\ {1} & {1}\end{bmatrix}, A_{3}=\begin{bmatrix}{1} & {1} & {0} \\ {1} & {1} & {1} \\ {0} & {1} & {1}\end{bmatrix}, A_{4}=\begin{bmatrix}{1} & {1} & {0} & {0} \\ {1} & {1} & {1} & {0} \\ {0} & {1} & {1} & {1} \\ {0} & {0} & {1} & {1}\end{bmatrix},\cdots, 寻找其行列式值的规律$$
显然$A_1,A_2$是容易求得的，$A_1=1,A_2=0$.

我们主要关注$A_3与A_4$，那么我们首先看$A_3$，利用代数余子式的方式，我们对第一行展开，会有:
$$显然直接按第一行展开有： A_3 =1\begin{vmatrix} 1&1\\1&1\end{vmatrix}-1\begin{vmatrix}1&1\\0& 1\end{vmatrix}=0-1=-1。$$我们把$\begin{vmatrix}1&1\\0& 1\end{vmatrix}按第一列展开，即\begin{vmatrix}1\end{vmatrix}=A_1$，那么我们会发现$A_3=A_2-A_1$
同样对于$A_4$，有:

$$A_4=\begin{vmatrix}1& 1& 0& 0\\1& 1& 1& 0\\0& 1& 1& 1\\0& 0& 1& 1\end{vmatrix}\stackrel{沿第一行展开}{=}1\begin{vmatrix}1& 1& 0\\1& 1& 1\\0& 1& 1\end{vmatrix}-1\begin{vmatrix}1& 1& 0\\0& 1& 1\\0& 1& 1\end{vmatrix}=A_3-1\times\begin{vmatrix} 1& 1\\ 1& 1\end{vmatrix}=A_3-A_2$$
那么我们根据数学归纳法，可以得到元素为1的上三角行列式值的递推式:
$$A_n=A_{n-1}-A_{n-2}$$
由此规律,易得 $|{A}_5|=0$, $|{A}_6|=1$, $|{A}_7|=1$, $|{A}_8|=0$, 到这里我们发现:由1组成的$n$阶三对角矩阵的行列式值从1阶开始按照1,0,-1,-1,0,1循环,周期为6。
