# Econometrics-using-Python
## 最小二乘估计 $\left(OLS\right)$

对于样本集 $\left\{\left(x_i,y_i\right)\right\}_{i=1,2,...,n}$ ，线性回归模型为：$y=X\beta+\varepsilon$ 。其中， $\beta=\left[\beta_1,\beta_2,...,\beta_K\right]^T$ ， $X=\left[x_1,x_2,\ldots,x_n\right]^T$ ， $y=\left[y_1,y_2,\ldots,y_n\right]^T$ 。
用 $b$ 估计 $\beta$ ，可得： $y=Xb+e$ 。
最优化问题为： ${min}\below b{e^Te}$ 
 $e^Te=\left(y-Xb\right)^T\left(y-Xb\right)=y^Ty-y^TXb-\left(Xb\right)^Ty+\left(Xb\right)^TXb=y^Ty-2y^TXb+b^TX^TXb$ 
y^TXb=\sum_{i=1}^{n}y_i{x_i}^Tb=\sum_{i=1}^{n}y_i\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)
\frac{\partial y^TXb}{\partial b}=\left[\frac{\partial y^TXb}{\partial b_1},\frac{\partial y^TXb}{\partial b_2},\ldots,\frac{\partial y^TXb}{\partial b_K}\right]^T=\left[\sum_{i=1}^{n}{y_ix_{i,1}},\sum_{i=1}^{n}{y_ix_{i,2}},\ldots,\sum_{i=1}^{n}{y_ix_{i,K}}\right]^T=\sum_{i=1}^{n}\left[y_ix_{i,1},y_ix_{i,2},..,y_ix_{i,K}\right]^T=\sum_{i=1}^{n}{y_i\left[x_{i,1},x_{i,2},..,x_{i,K}\right]^T}=\sum_{i=1}^{n}{y_i{x_i}^T}=y^TX
b^TX^TXb=b^T\sum_{i=1}^{n}{x_i{x_i}^T}b=\sum_{i=1}^{n}{b^Tx_i{x_i}^Tb}=\sum_{i=1}^{n}{\left({x_i}^Tb\right)^T{x_i}^Tb}=\sum_{i=1}^{n}\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)^2
\frac{\partial b^TX^TXb}{\partial b}=\left[\frac{\partial b^TX^TXb}{\partial b_1},\frac{\partial b^TX^TXb}{\partial b_2},\ldots,\frac{\partial b^TX^TXb}{\partial b_K}\right]^T=\left[2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,1}},2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,2}},\ldots,2\sum_{i=1}^{n}{\left(\sum_{k=1}^{K}{x_{i,k}b_k}\right)x_{i,K}}\right]^T=\left[2\sum_{i=1}^{n}{{x_i}^Tbx_{i,1}},2\sum_{i=1}^{n}{{x_i}^Tbx_{i,2}},\ldots,2\sum_{i=1}^{n}{{x_i}^Tbx_{i,K}}\right]^T=2\left[\sum_{i=1}^{n}{{x_i}^Tbx_{i,1}},\sum_{i=1}^{n}{{x_i}^Tbx_{i,2}},\ldots,\sum_{i=1}^{n}{{x_i}^Tbx_{i,K}}\right]^T=2\sum_{i=1}^{n}\left[{x_i}^Tbx_{i,1},{x_i}^Tbx_{i,2},\ldots,{x_i}^Tbx_{i,K}\right]^T=2\sum_{i=1}^{n}{{x_i}^Tb\left[x_{i,1},x_{i,2},\ldots,x_{i,K}\right]^T}=2\sum_{i=1}^{n}{{x_i}^Tbx_i}=2\sum_{i=1}^{n}{x_i{x_i}^Tb}=2\left(\sum_{i=1}^{n}{x_i{x_i}^T}\right)b=2X^TXb
\frac{\partial e^Te}{\partial b}=-2y^TX+2X^TXb=0
正规方程组：
X^TXb=y^TX
解得:
b^*=\left(X^TX\right)^{-1}y^TX
H_{e^Te}\left(b\right)=\frac{\partial2X^TXb}{\partial b}=2X^TX
\forall a\neq0，a^TX^TXa={(Xa)}^TXa>0，故有H_{e^Te}\left(b\right)矩阵正定，可得{argmin}\below b{e^Te}=b^*。
