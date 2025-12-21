# Econometrics-using-Python
## 一、最小二乘估计 $\left(OLS\right)$

对于样本集 $\\{\left(x_i,y_i\right)_{i=1,2,...,n}\\}$ ，线性回归模型为： $y=X\beta+\varepsilon$ 。其中， $𝛽=\left[𝛽_1,𝛽_2,...,𝛽_𝐾\right]^𝑇$ ， $𝑋=\left[𝑥_1,𝑥_2,…,𝑥_𝑛\right]^𝑇$ ， $𝑦=\left[𝑦_1,𝑦_2,…,𝑦_𝑛\right]^𝑇$ 。

用 $b$ 估计 $\beta$ ，可得： $y=\hat{y}+e=Xb+e$ 。

最优化问题为： $\min\limits_b{e^Te}$

$e^Te=\left(y-Xb\right)^T\left(y-Xb\right)=y^Ty-y^TXb-\left(Xb\right)^Ty+\left(Xb\right)^TXb=y^Ty-2y^TXb+b^TX^TXb$
 
$y^TXb=\sum_{i=1}^{n}y_i{x_i}^Tb=\sum_{i=1}^{n}y_i\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)$

$\frac{\partial y^TXb}{\partial b}=\left[\frac{\partial y^TXb}{\partial b_1},\frac{\partial y^TXb}{\partial b_2},\ldots,\frac{\partial y^TXb}{\partial b_K}\right]^T=\left[\sum_{i=1}^{n}{y_ix_{i,1}},\sum_{i=1}^{n}{y_ix_{i,2}},\ldots,\sum_{i=1}^{n}{y_ix_{i,K}}\right]^T=\sum_{i=1}^{n}\left[y_ix_{i,1},y_ix_{i,2},..,y_ix_{i,K}\right]^T=\sum_{i=1}^{n}{y_i\left[x_{i,1},x_{i,2},..,x_{i,K}\right]^T}=\sum_{i=1}^{n}{y_ix_i}=X^Ty$

$b^TX^TXb=b^T\sum_{i=1}^{n}{x_i{x_i}^T}b=\sum_{i=1}^{n}{b^Tx_i{x_i}^Tb}=\sum_{i=1}^{n}{\left({x_i}^Tb\right)^T{x_i}^Tb}=\sum_{i=1}^{n}\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)^2$

$\frac{\partial b^TX^TXb}{\partial b}=\left[\frac{\partial b^TX^TXb}{\partial b_1},\frac{\partial b^TX^TXb}{\partial b_2},\ldots,\frac{\partial b^TX^TXb}{\partial b_K}\right]^T=\left[2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,1}},2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,2}},\ldots,2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,K}}\right]^T=\left[2\sum_{i=1}^{n}{{x_i}^Tbx_{i,1}},2\sum_{i=1}^{n}{{x_i}^Tbx_{i,2}},\ldots,2\sum_{i=1}^{n}{{x_i}^Tbx_{i,K}}\right]^T=2\left[\sum_{i=1}^{n}{{x_i}^Tbx_{i,1}},\sum_{i=1}^{n}{{x_i}^Tbx_{i,2}},\ldots,\sum_{i=1}^{n}{{x_i}^Tbx_{i,K}}\right]^T=2\sum_{i=1}^{n}\left[{x_i}^Tbx_{i,1},{x_i}^Tbx_{i,2},\ldots,{x_i}^Tbx_{i,K}\right]^T=2\sum_{i=1}^{n}{{x_i}^Tb\left[x_{i,1},x_{i,2},\ldots,x_{i,K}\right]^T}=2\sum_{i=1}^{n}{{x_i}^Tbx_i}=2\sum_{i=1}^{n}{x_i{x_i}^Tb}=2\left(\sum_{i=1}^{n}{x_i{x_i}^T}\right)b=2X^TXb$

$\frac{\partial e^Te}{\partial b}=-2X^Ty+2X^TXb=0$

正规方程组： $X^TXb=X^Ty$

$OLS$ 假设 $r\left(X\right)=K$ ，可得： $r\left(X^TX\right)=K$

解得： $b^\ast=\left(X^TX\right)^{-1}X^Ty$

计算二阶海塞矩阵： $H_{e^Te}\left(b\right)=\frac{\partial2X^TXb}{\partial b}=2X^TX$

$\forall a\neq0$ ， $a^TX^TXa={(Xa)}^TXa>0$ ，故有 $H_{e^Te}\left(b\right)$ 正定，可得 $\mathop{argmin}\limits_b{e^Te}=b^*$ 。

## 二、 $\left(OLS\right)$ 性质
1.无偏性：

$OLS$ 假设 $E(\varepsilon|X)=0$

$E\left(\left(X^TX\right)^{-1}X^T\left(X\beta+\varepsilon\right)\middle| X\right)$

$E\left(b^\ast\middle| X\right)=E\left(\left(X^TX\right)^{-1}X^Ty\middle| X\right)=E\left(\left(X^TX\right)^{-1}X^T\left(X\beta+\varepsilon\right)\middle| X\right)=\beta+E\left(\left(X^TX\right)^{-1}X^T\varepsilon\middle| X\right)=\beta+\left(X^TX\right)^{-1}X^TE\left(\varepsilon\middle| X\right)=\beta$

2.一致性：

$OLS$ 假设 $E(\varepsilon|X)=0$

$^\ast=\left(X^TX\right)^{-1}X^Ty=\left(X^TX\right)^{-1}X^T\left(X\beta+\varepsilon\right)=\beta+\left(X^TX\right)^{-1}X^T\varepsilon$

$b^\ast-\beta=\left(X^TX\right)^{-1}X^T\varepsilon=\left(\frac{X^TX}{n}\right)^{-1}\left(\frac{X^T\varepsilon}{n}\right)$

$\left(\frac{X^TX}{n}\right)^{-1}=\left(\frac{1}{n}\sum_{i=1}^{n}{x_i{x_i}^T}\right)^{-1}$

假定 $x_1,x_2,\ldots,x_n$ 为抽样自同一总体的独立同分布样本， $E(x_i)=\mu，Var(x_i)=\sigma^2$ ，当 $n\rightarrow\infty$ 时， $\frac{1}{n}\sum_{i=1}^{n}x_i$ 依概率收敛于 $\mu$ 。

记 $E\left(x_i{x_i}^T\right)=Q$

$\mathop{Plim}\limits_{n\rightarrow\infty}{\frac{X^TX}{n}}=Q$

$\mathop{Plim}\limits_{n\rightarrow\infty}{\left(\frac{X^TX}{n}\right)^{-1}}=Q^{-1}$

$\frac{X^T\varepsilon}{n}=\frac{1}{n}\sum_{i=1}^{n}{x_i\varepsilon_i}$

$E\left(x_i\varepsilon_i\right)=E\left(E\left(x_i\varepsilon_i|X\right)\right)=E\left(x_iE\left(\varepsilon_i|X\right)\right)=0$

$\mathop{Plim}\limits_{n\rightarrow\infty}{\frac{X^T\varepsilon}{n}}=0$

$\mathop{Plim}\limits_{n\rightarrow\infty}\left(b^\ast-\beta\right)=0$

$\mathop{Plim}\limits_{n\rightarrow\infty}b^\ast=\beta$

3.有效性（最小方差性）：

$OLS$ 假设 $Var(\varepsilon|X)=\sigma^2I$

设 $\widetilde{b}$ 为线性回归 $y=X\beta+\varepsilon$ 的又一线性无偏估计量，即 $\widetilde{b}=Cy=\left(D+\left(X^TX\right)^{-1}X^T\right)y$

$E\left(\widetilde{b}\middle| X\right)=E\left(\left(D+\left(X^TX\right)^{-1}X^T\right)y\middle| X\right)=E\left(\left(D+\left(X^TX\right)^{-1}X^T\right)\left(X\beta+\varepsilon\right)\middle| X\right)=E\left(DX\beta+\beta+D\varepsilon+\left(X^TX\right)^{-1}X^T\varepsilon\middle| X\right)=\left(I+DX\right)\beta$

故： $I+DX=I$ ，即 $DX=0$

$Var\left(\widetilde{b}\middle| X\right)=Var\left(DX\beta+\beta+D\varepsilon+\left(X^TX\right)^{-1}X^T\varepsilon\middle| X\right)=DVar\left(\varepsilon\middle| X\right)D^T+\left(X^TX\right)^{-1}X^TVar\left(\varepsilon\middle| X\right)X\left(X^TX\right)^{-1}=\sigma^2DD^T+\sigma^2\left(X^TX\right)^{-1}X^TX\left(X^TX\right)^{-1}=\sigma^2DD^T+\sigma^2\left(X^TX\right)^{-1}$

$Var\left(b^\ast\middle| X\right)=Var\left(\left(X^TX\right)^{-1}X^T\varepsilon\middle| X\right)=\left(X^TX\right)^{-1}X^TVar\left(\varepsilon\middle| X\right)X\left(X^TX\right)^{-1}=\sigma^2\left(X^TX\right)^{-1}X^TX\left(X^TX\right)^{-1}=\sigma^2\left(X^TX\right)^{-1}$

因为 $DD^T$ 正定，故 $Var\left(b^\ast\middle| X\right)$ 为最小方差。

4.渐进正态性：

$b^\ast-\beta=\left(X^TX\right)^{-1}X^T\varepsilon=\left(\frac{X^TX}{n}\right)^{-1}\frac{X^T\varepsilon}{n}=\left(\frac{1}{n}\sum_{i=1}^{n}{x_i{x_i}^T}\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^{n}{x_i\varepsilon_i}\right)$

假定 $x_1,x_2,\ldots,x_n$ 为抽样自同一总体的独立同分布样本， $E\left(x_i\right)=\mu$ ， $Var\left(x_i\right)=\sigma^2$ ， 当 $n\rightarrow\infty$ 时， $\sqrt n\left(\frac{1}{n}\sum_{i=1}^{n}x_i-\mu\right)$ 依分布收敛于 $N\left(0,\sigma^2\right)$ 。

记 $E\left(x_i{x_i}^T\right)=Q$ ：

$Var\left(x_i\varepsilon_i\right)=E\left(\left(x_i\varepsilon_i-E\left(x_i\varepsilon_i\right)\right)\left(x_i\varepsilon_i-E\left(x_i\varepsilon_i\right)\right)^T\right)=E\left(x_i{\varepsilon_i}^2{x_i}^T\right)=E\left(x_iE\left({\varepsilon_i}^2\middle| X\right){x_i}^T\right)=\sigma^2E\left(x_i{x_i}^T\right)=\sigma^2Q$

$\sqrt n\left(\frac{1}{n}\sum_{i=1}^{n}{x_i\varepsilon_i}\right)\stackrel{d}{\rightarrow}N\left(0,\sigma^2Q\right)$

$\mathop{Plim}\limits_{n\rightarrow\infty}{\left(\frac{1}{n}\sum_{i=1}^{n}{x_i{x_i}^T}\right)^{-1}}=Q^{-1}$

$\sqrt n\left(\frac{1}{n}\sum_{i=1}^{n}{x_i{x_i}^T}\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^{n}{x_i\varepsilon_i}\right)\stackrel{d}{\rightarrow}N\left(0,\sigma^2\left(Q^{-1}\right)^TQQ^{-1}\right)$

$\left(Q^{-1}\right)^T=ExixiT-1T=ExixiTT-1=ExixiT-1=Q-1$

$\sqrt n\left(\frac{1}{n}\sum_{i=1}^{n}{x_i{x_i}^T}\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^{n}{x_i\varepsilon_i}\right)\stackrel{d}{\rightarrow}N\left(0,\sigma^2Q^{-1}\right)$

$\sqrt n\left(b^\ast-\beta\right)\stackrel{d}{\rightarrow}N\left(0,\sigma^2Q^{-1}\right)$
