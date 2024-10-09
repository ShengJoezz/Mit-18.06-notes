---
author: joe
rating: 8
content: 介绍了行列式最基本的几个三个运算法则和基于此的几个二级的复杂运算法则
coverage: p18
tags:
  - linear_algebra
  - math
  - class
---

迄今为止,我们已经学习了很多关于长方矩阵的知识,现在,把注意力转向方阵,探讨两个大的话题:行列式和特征值,我们需要行列式的重要原因是求特征值。

行列式是跟每个方阵都有关的数,每个方阵都有与其相关的行列式,我们一般记为 $\det {A}$,或者写作 $|{A}|$。

行列式最早是应用在用来判断方程组是否有解,在矩阵被发明后,行列式就拥有了更多的性质和应用。行列式是一个神奇的数,一个数很难告诉你整个矩阵是什么样子的,但行列式把矩阵的尽可能多的信息就包含在其中了。就比如,矩阵可逆等价于行列式非零,行列式为零时矩阵是[[奇异矩阵]]的。从另外一个角度来理解,行列式从某种程度上代表了这个矩阵的特征,这是学习特征分解的前置概念。
# Section1.三个基本运算法则
##  性质 1: 对于单位矩阵 ${I}$,有 $\det {I} = 1$。
## 性质 2: 交换行,行列式的值的符号会相反。

对于性质二，我们考虑一个问题，是否会出现我换7行的结果与换10行的结果行列式的长相一样呢？
事实上，对于一个矩阵,我通过 $7$ 次换行得到一种置换,那么同样可以通过 $21$ 或 $23$ 次换行得到甚至是 $101$ 次的,不管具体是多少但这个置换一定是奇数次的而不会是偶数次的。
对于一个给定的矩阵,如果该矩阵是可逆的,那么该矩阵 $7$ 次行交换的结果,绝对不会与 $10$ 次行交换的结果相同。Gilbert教授没有给出具体的证明，这种置换的奇偶区分意味着行列式性质二是严谨且正确的,同时意味着行列式可以被严格定义。
```ad-note
title: 推断
根据性质1与2，我们可以知道[[Chapter03矩阵的转置，置换#Section1. 置换矩阵|置换矩阵]]的行列式是$1$或$-1$，当然具体是$1$还是$-1$我们需要考虑行交换的次数:
$$\det P = \begin{cases} 1 & \text{even (偶数)}\\ -1 & \text{odd (奇数)}  
\end{cases}$$

```
## 性质3:行拆分的性质
### 性质 3.a: 行列式按行提取出矩阵中的系数,也即 $\begin{vmatrix}ta&tb\\c&d\end{vmatrix}=t\begin{vmatrix}a&b\\c&d\end{vmatrix}$
### 性质 3.b: 行列式是一个线性函数,但这个线性单独反映在每一行上,也即: $\begin{vmatrix}a+k_1&b+k_2\\c&d\end{vmatrix}=\begin{vmatrix}a&b\\c&d\end{vmatrix}+\begin{vmatrix}k_1&k_2\\c&d\end{vmatrix}$,只要有一行是相同的，就可以对第一行进行任意的拆分组合。

```ad-note
title: 注意
我们注意性质3.b，这个线性特性针对的是单行而不是指整个矩阵，说$det{A+B}=det{A}+det{B}$是错误的。
$$ \begin{vmatrix}a+a'&b+b'&c+c'\\d&e&f\\g+g'&h+h'&i+i'\end{vmatrix}=\begin{vmatrix}a&b&c\\d&e&f\\g+g'&h+h'&i+i'\end{vmatrix}+\begin{vmatrix}a'&b'&c'\\d&e&f\\g+g'&h+h'&i+i'\end{vmatrix}=\begin{vmatrix}a&b&c\\d&e&f\\g&h&i\end{vmatrix}+\begin{vmatrix}a&b&c\\d&e&f\\g'&h'&i'\end{vmatrix}+\begin{vmatrix}a'&b'&c'\\d&e&f\\g&h&i\end{vmatrix}+\begin{vmatrix}a'&b'&c'\\d&e&f\\g'&h'&i'\end{vmatrix} $$
上面这个实例说明，针对行是能够进行任意的拆分的，如果要拆分两行分别成为两部分，就是2*2=4个部分，这也是好理解的。

```
```ad-note
title: 关于行列式线性性质的行拆分
我们知道$\begin{vmatrix}a+k_1&b+k_2\\c&d\end{vmatrix}=\begin{vmatrix}a&b\\c&d\end{vmatrix}+\begin{vmatrix}k_1&k_2\\c&d\end{vmatrix}$
这是只有单行参与的行拆分，实际上我们可以多行同时参与,例如$$\begin{vmatrix}a+k_1&b+k_2\\c+k_3&d+k_4\end{vmatrix}=\begin{vmatrix}a&b\\c&d\end{vmatrix}+\begin{vmatrix}a&b\\k_3&k_4\end{vmatrix}+\begin{vmatrix}k_1&k_2\\c&d\end{vmatrix}+\begin{vmatrix}k_1&k_2\\k_3&k_4\end{vmatrix}$$
也就是说我们把每一行确定一个组合.第一行:$\begin{vmatrix}a&b\end{vmatrix}与\begin{vmatrix}k_1&k_2\end{vmatrix}$
第二行:$\begin{vmatrix}c&d\end{vmatrix}与\begin{vmatrix}k_3&k_4\end{vmatrix}$，接着排列组合就好。

```
# Section2.推导的七个性质

## 性质4.如果两行相等，行列式等于0

证明:
根据[[#性质 2 交换行,行列式的值的符号会相反。]]，我们考虑交换过相同两行的行列式$detA'$，由于交换的是相同的两行，$detA与detA'完全一致$，则有$detA=detA'=-detA$，因此$detA=0$.

## 性质5.从 $k$ 行中减去第 $i$ 行的 $l$ 倍,行列式不变。

这条性质非常有用，它告诉我们，一个行列式未知的矩阵，消元得到其化简的形式，比如说上三角形式，**行列式并不因消元而改变**。

证明:
根据[[#性质 3.b 行列式是一个线性函数,但这个线性单独反映在每一行上,也即 $ begin{vmatrix}a+k_1&b+k_2 c&d end{vmatrix}= begin{vmatrix}a&b c&d end{vmatrix}+ begin{vmatrix}k_1&k_2 c&d end{vmatrix}$,只要有一行是相同的，就可以对第一行进行任意的拆分组合。|性质3.b与3.a]]，我们可以把减去的行做拆分，用简单的三维矩阵举例，即:
$$\begin{vmatrix}a+kg&b+kh&c+ki\\d&e&f\\g&h&i\end{vmatrix}=\begin{vmatrix}a&b&c\\d&e&f\\g&h&i\end{vmatrix}+k\begin{vmatrix}g&h&i\\d&e&f\\g&h&i\end{vmatrix}$$
而根据[[#性质4.如果两行相等，行列式等于0]]，可知$k\begin{vmatrix}g&h&i\\d&e&f\\g&h&i\end{vmatrix}=0$,即证。
## 性质6.如果方阵的某一行全为$0$，那么其行列式值也为$0$。

证明:
根据[[#性质 3.a 行列式按行提取出矩阵中的系数,也即 $ begin{vmatrix}ta&tb c&d end{vmatrix}=t begin{vmatrix}a&b c&d end{vmatrix}$|性质3.a]],可以知道$$\begin{vmatrix}0&0\\c&d\end{vmatrix}=0\begin{vmatrix}a&b\\c&d\end{vmatrix}=0$$
## 性质7.上三角矩阵对应的行列式的值等于其对角线上元素的乘积。


证明:
设有上三角矩阵 $U=\begin{vmatrix} d_1 & * & \cdots & *\\ 0 & d_2 & \cdots & *\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & d_n \end{vmatrix}$,则 $\det U = d_1d_2\cdots d_n$。
使用[[#性质5.从 $k$ 行中减去第 $i$ 行的 $l$ 倍,行列式不变。|性质5(某行减去其余行行列式不变)]],从最后一行开始,将对角元素上方的 $*$ 元素依次变为0,则可以得到型为$D=\begin{vmatrix} d_1 & 0 & \cdots & 0\\ 0 & d_2 & \cdots & 0\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & d_n \end{vmatrix}$的对角行列式,
再使用[[#性质 3.a 行列式按行提取出矩阵中的系数,也即 $ begin{vmatrix}ta&tb c&d end{vmatrix}=t begin{vmatrix}a&b c&d end{vmatrix}$|行提出t可以放在整个行列式外]]将对角元素提出得到$d_n\cdots d_1\begin{vmatrix} 1 & 0 & \cdots & 0\\ 0 & 1 & \cdots & 0\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & 1 \end{vmatrix}$,再结合[[#性质 1 对于单位矩阵 ${I}$,有 $ det {I} = 1$。|detI=1]],得证。
```ad-note
需要补充说明的是,从矩阵$A$化简到$U$可能中间除了消元还有换行过程,故矩阵$A$的真正的行列式可能不与$U$的行列式同号,这是由[[#性质 2 交换行,行列式的值的符号会相反。]]决定的。

```

## 性质 8: 当且仅当 ${A}$ 是奇异矩阵(不可逆)时,$\det {A}=0$。

$$奇异矩阵(不好的矩阵所以不可逆)\rightarrow不可逆\rightarrow秩小于n\rightarrow行秩小于n,主元不满，有0行\rightarrow根据性质6可知detA=0$$
```ad-note
title: 一个证明
我们知道 $\begin{vmatrix} a & c\\ b & d \end{vmatrix}=ad-bc$
我们可以验证一下,考虑二阶矩阵的情况, $\begin{vmatrix} a & c\\ b & d \end{vmatrix}$ $\rightarrow$ $\begin{vmatrix} a & 0\\ b & d- \frac{c}{a}b  \end{vmatrix} = ad - bc$

此时求行列式的值,发现对角元素确实等于主元相乘。

```
## 性质9:$detAB=detAdetB$

```ad-note
title: 有趣的点
行列式不具备线性性质，不具备加法性质，但却具备乘法性质。

```
使用这一性质,我们就很好求得
* ${A}$ 的逆矩阵的行列式: $\det {I} = \det {A}^{-1}{A} = \det {A}^{-1}\det {A}=1$ $\therefore \det {A}^{-1} = \frac{1}{\det {A}}$
*  $\det {A}^2 = (\det {A})^2$:矩阵平方的行列式等于矩阵行列式的平方。
* $\det k{A} = detkI\times detA=k^n\det {A}$:$\ k$ 为常数, ${A}$ 为 $n$ 阶的矩阵,这里相当于提出了每一行的 $k$,这个式子就像是在求体积一样，改变立方体的边长为原来的$k$倍，但是体积变为原来的$k^3$倍。

给出一个粗略的证明，任何的矩阵都可以化为上三角矩阵，并且根据[[#性质5.从 $k$ 行中减去第 $i$ 行的 $l$ 倍,行列式不变。|性质5]]，可以知道行列式不变，那么我们考虑:
$$AB=\begin{bmatrix} a_1 & 0 & \cdots & 0\\ 0 & a_2 & \cdots & 0\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & a_n \end{bmatrix} \begin{bmatrix} b_1 & 0 & \cdots & 0\\ 0 & b_2 & \cdots & 0\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & b_n \end{bmatrix}=
\begin{bmatrix} a_1b_1 & 0 & \cdots & 0\\ 0 & a_2b_2 & \cdots & 0\\ \vdots & \vdots & \ddots & \vdots\\ 0 & 0 & \cdots & a_nb_n  
\end{bmatrix}$$
而对角矩阵其行列式就等于对角元素之乘积,所以 $\det(AB) = (a_1b_1)\cdots(a_nb_n) = (\prod_{i=1}^n a_i)(\prod_{i=1}^n b_i) = \det(A)\det(B)$。(这里我们注意到，哪怕$a_i/b_i=0$依然成立)

## 性质 10: $\det {A}^T = \det {A}$.

证明:将矩阵化为${LU}$的形式有:${A}^T={U}^T{L}^T$, ${A}={LU}$。那么由[[#性质9 $detAB=detAdetB$|行列式的乘法法则]]可得 $|{A}^T|=|{U}^T{L}^T|=|{U}^T||{L}^T|$, $|{A}|=|{LU}|=|{L}||{U}|$,又因为${L},{U}$的行列式并不因为转置而改变(因为$L$是上三角矩阵,$U$是下三角矩阵),于是得证。