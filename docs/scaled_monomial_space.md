# 缩放单项式空间
<font size = 1>

|次幂 |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |
|--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |
|$k=0$: |    |   |   |    |    | 0: $1$ |    |    |    |    |    |
|$k=1$: |     |    |    |    |  1:$\bar x$ |    |  2:$\bar y$
|$k=2$: |     |    |    |  3:${\bar x}^2$ |    |  4:$\bar x\bar y$ |    | 5:${\bar y}^2$
|$k=3$: |     |    |  6:${\bar x}^3$ |    |  7:${\bar x}^2\bar y$ |    | 8:$\bar x{\bar y}^2$ |    |  9:${\bar y}^3$
|$k=4$: |     |  10:${\bar x}^4$ |    |  11:${\bar x}^3\bar y$ |    | 12:${\bar x}^2{\bar y}^2$ |    |  13:${\bar x}{\bar y}^3$ |    | 14:${\bar y}^4$
|$k=5$: |   15:${\bar x}^5$ |    |  16:${\bar x}^4\bar y$ |    | 17:${\bar x}^3{\bar y}^2$ |    | 18:${\bar x}^2{\bar y}^3$ |    | 19:${\bar x}{\bar y}^4$ |    |  20:${\bar y}^5$|
</font>

令 $K$ 是一个 $R^2$ 上的多边形, 面积为 $|K|$, 尺寸为 $h_K = \sqrt{|K|}$, 
重心为 $\boldsymbol x_K = (x_K, y_K)$. 定义:

$$
\bar x = \frac{x - x_K}{h_K}, \quad \bar y = \frac{y - y_K}{h_K}
$$

记 $\boldsymbol \alpha = (\alpha_0, \alpha_1)$ 为任一二重非负整数指标, 
则 $K$ 上的 **缩放单项式** 可表示为:

$$
\boldsymbol m_{\boldsymbol \alpha} = \bar{x}^{\alpha_0} \bar{y}^{\alpha_1}
$$

如上表. 记 $m_k$ 为 $k$ 次多项式组成的向量函数, 按上表中方式排序.

$$
\boldsymbol m_k = [m_0, m_1, ..., m_{n_k-1}]
$$

则显然有:

$$
\boldsymbol m_k[:-2] = \bar x * \boldsymbol m_{k-1}, \quad 
\boldsymbol m_k[-1] = \bar y * \boldsymbol m_{k-1}[-1]
\tag{1}
$$

$K$ 上的 $p$ 次多项式空间的基函数即为: 
$m_0 \cup m_1 \cup ...\cup m_p$, 对于使用 $(1)$ 式可计算每个基函数.

基函数导数计算公式:

$$
\frac{\partial \boldsymbol m_k}{\partial x}[:-2] = 
[k, k-1, ..., 1]*\boldsymbol m_{k-1}, \quad 
\frac{\partial \boldsymbol m_k}{\partial x}[-1] = 0
$$

$$
\frac{\partial \boldsymbol m_k}{\partial y}[1:] = 
[1, 2, ..., k]*\boldsymbol m_{k-1}, \quad
\frac{\partial \boldsymbol m_k}{\partial x}[0] = 0
$$

二阶导计算公式:

$$
\frac{\partial^2 \boldsymbol m_k}{\partial x^2}[:-3] = 
[k*(k-1), (k-1)(k-2), ..., 2*1]*\boldsymbol m_{k-1}, 
$$

$$
\frac{\partial^2 \boldsymbol m_k}{\partial x^2}[-2:] = 0
$$

$$
\frac{\partial^2 \boldsymbol m_k}{\partial y^2}[2:] = 
[1*2, 2*3, ..., (k-1)*k]*\boldsymbol m_{k-1},
$$

$$
\frac{\partial^2 \boldsymbol m_k}{\partial y^2}[:2] = 0
$$

$$
\frac{\partial^2 \boldsymbol m_k}{\partial x \partial y}[1:-2] = 
[1*(k-1), 2*(k-2), ..., (k-1)*1]*\boldsymbol m_{k-1},
$$

$$
\frac{\partial^2 \boldsymbol m_k}{\partial x \partial y}[[0, -1]] = 0
$$



























