---
author: joe
rating: 8
content: 解线性方程组的视角转化为向量的线性组合，从而利用矩阵语言的Gauss消元法行变换来消元得解
coverage: p1-p3
tags:
  - math
  - linear_algebra
  - class
---

# Section1.引入 (线性方程组的求解)
## Subsection1.线性方程组视角转换

我们在考虑线性方程组的求解中，往往是以`row picture`(行图像)的视角来讲述的，例如我们针对如下的方程组:

$$
\left\{
\begin{aligned}
2x - y &= 0 \\
-x + 2y &= 3
\end{aligned}
\right.
$$

视为AX=b

我们可以将系数矩阵A写为:

$$
\begin{bmatrix}
2 & -1\\
-1 & 2\\
\end{bmatrix}
$$

我们以`row picture`的视角来理解，实际上就是在考虑两个方程在平面直角坐标系中的相交点。

实际上，我们可以用另一个视角。即`column picuture`的视角来理解,把系数矩阵考虑为两个列向量:

$$
\alpha=
\begin{bmatrix}
2 \\
-1\\
\end{bmatrix}
$$
$$
\beta=
\begin{bmatrix}
-1\\
2\\
\end{bmatrix}
$$
那么我们就考虑:
$$x\alpha+y\beta=b$$
也就是把b视为α与β的线性组合(`linear combination`),这是一种视角的转变 ^92b2a0

## Subsection2.矩阵都对任意b的有解性

答案是否定的，这对系数矩阵`A`的性质有一定要求，即关注`A`是否是[[奇异矩阵]]或者叫做是否是[[可逆矩阵]]，事实上，这个与[[Chapter02矩阵的求逆与A=LU分解法]]的可行性也是有关系的

## Subsecion3.消元法

我们需要确定主元，n个变量的方程组当然是n个主元，例如我们这里拿三维的举例，有:

$$
\left\{
\begin{aligned}
x +2y+z &= 2 \\
3x + 8y+z &= 12\\
4y+z=2\\
\end{aligned}
\right.
$$

系数矩阵就可以写为:

$$A=
\begin{bmatrix}
1& 2 & 1\\
3 & 8& 1\\
0& 4 & 1\\
\end{bmatrix}
$$
我们确定三个主元并进行 行之间的相减，可以得到:
$$U=
\begin{bmatrix}
1& 2 & 1\\
0 & 2& -2\\
0& 0 & 5\\
\end{bmatrix}
$$

^5816e6

这样我们就把主元之前的系数全部写为0了。现在的工作是回代，即`augument matrix`(增广矩阵),将方程右侧的三个数字进行行向量相同的变化。

同样的，不是所有的矩阵都能使用消元法得到解，我们要求**主元是不能为0**的，这里如果最后一个方程是`4y-4z=2`，这个方程组仍旧无法得到我们想要的结果。这个就是我们上面提到的[[奇异矩阵]]与[[可逆矩阵]]。

## Subsecion4.利用初等矩阵进行的Gauss消元(行变换)

在[[Chapter01线性方程组的求解#^92b2a0|解向量是向量的线性组合]]中我们已经窥见一点，即:

$$
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix} \cdot \begin{bmatrix}
a\\
b\\
\end{bmatrix} = a\cdot \begin{bmatrix}
1\\3
\end{bmatrix}+b\cdot\begin{bmatrix}
2\\4
\end{bmatrix}
$$
在某个矩阵进行右乘的时候，实际上是对列向量进行变换。

同样的，我们把这个性质考虑到对行向量上，那么我们可以得到:

$$
\begin{bmatrix}
a &b\\
\end{bmatrix} \cdot 
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}  = a\cdot \begin{bmatrix}
1&2
\end{bmatrix}+b\cdot\begin{bmatrix}
3&4
\end{bmatrix}
$$

我们将他稍微做一点点拓展，我们这里的结果只有一个向量，我们把左乘的矩阵做一点拓展，就可以得到一个类似的结果:

$$
\begin{bmatrix}
a &b\\
c&d\\
\end{bmatrix} \cdot 
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}  = 
\begin{bmatrix}
a\cdot\begin{bmatrix}
1&2
\end{bmatrix}+b\cdot\begin{bmatrix}
3&4
\end{bmatrix}\\
c\cdot\begin{bmatrix}
1&2
\end{bmatrix}+d\cdot\begin{bmatrix}
3&4
\end{bmatrix}
\end{bmatrix}
$$

而我们注意到，我们上面在做[[Chapter01线性方程组的求解#Subsecion3.消元法|Guass消元]]的时候，对于系数矩阵做的正是行变化，而这些行变化我们也正是我们在探究的，也就是说，我们可以用**矩阵的方式去书写消元的过程**

例如，我们在上面这个例子中，我们进行的步骤是，先消去第二行第一列的内容，我们做的动作是，第一行第三行保持不动，第二行减去第一行乘以3，我们用矩阵的形式来描述就是:

$$\begin{bmatrix}
1&0&0\\
-3&1&0\\
0&0&1\\
\end{bmatrix}\times A$$

这个矩阵我们称之为[[初等矩阵]]，由于是消除的第二行第一列的元素，我们写作$E_{21}$ ,同样的,我们可以做$E_{32}$ ，这样就可以形成最后的矩阵$U$ .即:

$$E_{32}E_{31}E_{21}A=U$$

如果我们只想用一个矩阵来说明这个变化过程，只需要把$E_{32}E_{31}E_{21}$ 计算出来即可.这说明矩阵乘法满足结合律.

## Subsection5.矩阵乘法的几个视角

* 最简单的运算法则，即所求的最终矩阵的某个元素$a_{ij}$ 是$\sum_{k=1}^na_{ik}b_{kj}$ ;

* [[Chapter01线性方程组的求解#Subsecion4.利用初等矩阵进行的Gauss消元(行变换)|行变化]]的视角，即
$$
\begin{bmatrix}
a &b\\
c&d\\
\end{bmatrix} \cdot 
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}  = 
\begin{bmatrix}
a\cdot\begin{bmatrix}
1&2
\end{bmatrix}+b\cdot\begin{bmatrix}
3&4
\end{bmatrix}\\
c\cdot\begin{bmatrix}
1&2
\end{bmatrix}+d\cdot\begin{bmatrix}
3&4
\end{bmatrix}
\end{bmatrix}
$$
* 列变化的视角，即
$$
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix} \cdot \begin{bmatrix}
a\\
b\\
\end{bmatrix} = a\cdot \begin{bmatrix}
1\\3
\end{bmatrix}+b\cdot\begin{bmatrix}
2\\4
\end{bmatrix}
$$
 要注意一点，行变化的视角我们一行一行地看，列变化视角我们一列一列地看.

* A每列×B每行的视角，

注意到column of A是$m\times1$ 的，而 row of B 是$1\times n$ 的，相乘结果是$m\times n$ 的，因此把每行每列拆开计算得到多个矩阵相加也可以。 ^5d559f

* 分块矩阵的视角

$$\begin{bmatrix}
A_{1} &A_{2}\\
A_{3} &A_{4}
\end{bmatrix}\times
\begin{bmatrix}
B_{1}&B_{2}\\
B_{3}&B_{4}
\end{bmatrix}$$
这个运算仍然符合最基础的矩阵乘法规则，例如最后结果的矩阵的左上角的结果应该为$A_{1}B_{1}+A_{2}B_{3}$ .

## 一个重要的视角:

我们理解矩阵比较容易接受的一个角度是:
列乘以行的视角:
$$Q^TQ=\begin{bmatrix} q_1 & q_2 & \cdots & q_k \end{bmatrix}\begin{bmatrix}q_1^T\\ q_2^T\\ \vdots\\ q_k^T   
\end{bmatrix}=\begin{bmatrix}q_1q_1^T+q_2q_2^T+\cdots+q_kq_k^T\end{bmatrix}$$即$$n\times k与k\times n=[n\times1与1\times n]=n\times n$$
还有行乘以列的视角:
$$Q^TQ=\begin{bmatrix} q_1^T\\ q_2^T\\ \vdots\\ q_n^T \end{bmatrix}\begin{bmatrix} q_1 & q_2 & \cdots & q_n  
\end{bmatrix}=\begin{bmatrix}q_1^Tq_1&q_1^Tq_2&\cdots&q_1^Tq_n\\&&\cdots\\&&\cdots\\q_n^Tq_1&q_n^Tq_2&\cdots&q_n^Tq_n\end{bmatrix}$$
即
$$n\times k与n\times k=[n\times n个1\times k与k\times 1]=n\times n$$
这两种都是正确的，实际上就是根据分块矩阵运算得来的。
==只要分块的矩阵符合运算规律，即分的块的前矩阵的列等于后面矩阵的行，就可以计算==