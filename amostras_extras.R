


```{r}
set.seed(13031998)

n=500

data_mean=data.frame(
  var0=1, # Intercepto média
  var1=rgamma(n,5,5/2),
  var2=rnorm(n,3,2),
  var3=0, # Intercepto variância
  var4=0
) %>% as.matrix

data_var=data.frame(
  var0=0, # Intercepto média
  var1=0,
  var2=0,
  var3=rep(1,n), # Intercepto variância
  var4=rpois(n,20)
) %>% as.matrix

beta_true=c(1,0.5,0.25,1,-0.04)
k=length(beta_true)

lin_pred_mean=as.matrix(data_mean)%*%beta_true
lin_pred_var=as.matrix(data_var)%*%beta_true

y=rnorm(n,lin_pred_mean,exp(-lin_pred_var/2))

# index=5
# plot(data_mean[,index]+data_var[,index],y)
```


```{r}
beta_mean=beta_mean_prior=rep(0,k)
beta_var=beta_var_prior=1*diag(k)

for(i in 1:n){
  FF=cbind(data_mean[i,],data_var[i,])
  x_mean=t(FF)%*%beta_mean
  Sigma=t(FF)%*%beta_var%*%FF
  mu1=x_mean[1]
  mu2=x_mean[2]
  s1=Sigma[1,1]
  s2=Sigma[2,2]
  rho=Sigma[1,2]/sqrt(s1*s2)
  val=rho*sqrt(s1/s2)

  x_var_adapt=(1-rho**2)*s1
  prec_prior=1/x_var_adapt
  f=function(x){
    x_mean_adapt=mu1+val*(log(x)-mu2)

    prec_obs=x
    x_var_update=1/(prec_obs+prec_prior)
    w1=1/(1+prec_prior/prec_obs)
    w2=1/(prec_obs/prec_prior+1)
    x_mean_update=w1*y[i]+w2*x_mean_adapt
    var_obs=1/prec_obs+x_var_adapt

    vals=dnorm(y[i],x_mean_adapt,sqrt(var_obs))*dlnorm(x,mu2,sqrt(s2))
    rbind(vals,
          log(x)*vals,
          (log(x)**2)*vals,
          x_mean_update*vals,
          (x_mean_update**2+x_var_update)*vals,
          (x_mean_update*log(x))*vals
    )

  }
  int_prev=cubintegrate(f,0,Inf,nVec=200,fDim=6)
  vals_prev=int_prev$integral
  c_prev=vals_prev[1]

  y_mean_prev=vals_prev[2]/c_prev
  y2_mean_prev=vals_prev[3]/c_prev
  y_s2_prev=y2_mean_prev-y_mean_prev**2

  x_mean_prev=vals_prev[4]/c_prev
  x2_mean_prev=vals_prev[5]/c_prev
  x_s2_prev=x2_mean_prev-x_mean_prev**2

  xy_mean_prev=vals_prev[6]/c_prev
  cov_prev=xy_mean_prev-x_mean_prev*y_mean_prev

  # sample=mvtnorm::rmvnorm(1000000,x_mean,Sigma)
  # mean((sample[,1]*sample[,2])*dnorm(y[i],sample[,1],exp(-sample[,2]/2)))/mean(dnorm(y[i],sample[,1],exp(-sample[,2]/2)))

  x_mean_new=c(x_mean_prev,y_mean_prev)
  Sigma_new=matrix(c(x_s2_prev,cov_prev,cov_prev,y_s2_prev),2,2)

  At <- beta_var %*% FF %*% ginv(Sigma)
  beta_mean <- beta_mean + At %*% (x_mean_new - x_mean)
  beta_var <- beta_var + At %*% (Sigma_new - Sigma) %*% t(At)
}
```


```{r}
log.like=function(beta){
  vals_mean=data_mean%*%beta
  vals_var=data_var%*%beta
  sum(dnorm(y,vals_mean,exp(-vals_var/2),log=TRUE))+dmvnorm(t(beta),beta_mean_prior,beta_var_prior,log=TRUE)
}

d_log_like=function(beta){calculus::derivative(function(beta){log.like(t(beta))},var=t(beta))}

mode=multiroot(d_log_like,start=beta_mean)$root

A=solve(-calculus::hessian(function(beta){log.like(t(beta))},var=t(mode)))
```

```{r}
beta_mean_ref=mode
beta_var_ref=A

N=5000
beta_sample=matrix(NA,k,N)
acep_list=matrix(0,k,N)
acep_p_list=matrix(0,k,N)
beta_unit=beta_mean_ref
cur_log_like=log.like(beta_unit)
beta_var_list=rep(NA,k)
for(j in 1:k){
  beta_var_list[j]=beta_var_ref[j,j]-t(beta_var_ref[j,-j])%*%ginv(beta_var_ref[-j,-j])%*%beta_var_ref[j,-j]
}
for(i in 1:N){
  for(j in 1:k){
    beta_mean_j=beta_mean_ref[j]+beta_var_ref[j,-j]%*%ginv(beta_var_ref[-j,-j])%*%(beta_unit[-j]-beta_mean_ref[-j])
    beta_var_j=beta_var_list[j]
    prop_j=rnorm(1,beta_mean_j,sqrt(beta_var_j))
    prop=beta_unit
    prop[j]=prop_j
    prop_log_like=log.like(prop)
    prop_log_val=dnorm(prop_j,beta_mean_j,sqrt(beta_var_j),log=TRUE)
    cur_log_val=dnorm(beta_unit[j],beta_mean_j,sqrt(beta_var_j),log=TRUE)
    acep=exp(prop_log_like-cur_log_like+cur_log_val-prop_log_val)
    acep_list[j,i]=acep
    if(runif(1)<acep){
      beta_unit[j]=prop_j
      cur_log_like=prop_log_like
      acep_list[j,i]=1
    }
  }
  beta_sample[,i]=beta_unit
}
beta_mean_true=rowMeans(beta_sample)
beta_var_true=var(t(beta_sample))
# mean(acep_list)
```



```{r}
data_plot=data.frame()
data_sample=pivot_longer(cbind(1:N,as.data.frame(t(beta_sample))),-1)
names(data_sample)=c('samp','index','value')
data_sample$index=substr(data_sample$index,2,2)
for(index in 1:k){
  x=seq(qnorm(1/N,beta_mean_true[index],sqrt(beta_var_true[index,index])),
        qnorm(1-1/N,beta_mean_true[index],sqrt(beta_var_true[index,index])),
        l=1000)
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,beta_mean[index],sqrt(beta_var[index,index])),
    label='Aprox. KL alt.',
    index=as.character(index)
  ))
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,mode[index],sqrt(A[index,index])),
    label='Aprox. Laplace',
    index=as.character(index)
  ))
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,beta_mean_true[index],sqrt(beta_var_true[index,index])),
    label='Aprox. KL orig.',
    index=as.character(index)
  ))

}

data_plot$index=as.factor(data_plot$index)
data_sample$index=as.factor(data_sample$index)
levels(data_sample$index)=levels(data_plot$index)=c('beta_0',
                                                    'beta_1',
                                                    'beta_2',
                                                    'beta_3',
                                                    'beta_4')
levels(data_sample$index)=levels(data_plot$index)=c('beta_0'=TeX('\\beta_0'),
                                                    'beta_1'=TeX('\\beta_1'),
                                                    'beta_2'=TeX('\\beta_2'),
                                                    'beta_3'=TeX('\\beta_3'),
                                                    'beta_4'=TeX('\\beta_4'))

(ggplot(data_plot %>% filter(label!='Aprox. Laplace'))+
    geom_histogram(aes(x=value,y=..density..),color='black',fill='#aaaaaa',bins=30,data=data_sample)+
    geom_line(aes(x=x,y=fx,color=label,linetype=label=='Aprox. KL orig.'))+
    scale_color_manual('',values=c('blue','red','green'))+
    scale_y_continuous(expand=c(0,0,0,0.5),limits=c(0,NA))+
    facet_wrap(.~index,scales='free')+
    guides(linetype='none')+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}

vals_mean_mean=data_mean%*%beta_mean
vals_var_mean=data_var%*%beta_mean

vals_mean_var=diag(data_mean%*%beta_var%*%t(data_mean))
vals_var_var=diag(data_var%*%beta_var%*%t(data_var))

vals_cov=diag(data_var%*%beta_var%*%t(data_mean))

mu0 <- vals_mean_mean + vals_cov
c0 <- exp(-vals_var_mean - vals_var_var / 2) / (vals_mean_var)
helper <- -3 + 3 * sqrt(1 + 2 * vals_var_var / 3)
# helper=Qt[2,2]
# print(c0)
# print(ft)
# print(Qt)
alpha <- 1 / helper
beta <- alpha * exp(-vals_var_mean - vals_var_var / 2)

nu <- 2 * alpha
sigma2 <- (beta / alpha) * (1 + 1 / c0)

y_pred_mean=mu0
y_pred_icl=qt(0.025,nu) * sqrt(sigma2) + mu0
y_pred_icu=qt(0.975,nu) * sqrt(sigma2) + mu0

vals_mean_sample=data_mean%*%beta_sample
vals_var_sample=data_var%*%beta_sample

y_sample=(rnorm(prod(dim(vals_mean_sample))) %>% matrix(n,N))*exp(-vals_var_sample/2)+vals_mean_sample

y_pred_mean_true=rowMeans(y_sample)
y_pred_icl_true=apply(y_sample,1,function(x){quantile(x,0.025)})
y_pred_icu_true=apply(y_sample,1,function(x){quantile(x,0.975)})


(ggplot()+
    geom_point(aes(x=1:n,y=y,color='Observações',fill='Observações'))+
    geom_point(aes(x=1:n,y=y_pred_mean,color='Aprox. KL',fill='Aprox. KL'),linetype='dashed',shape=3)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl,
                      ymax=y_pred_icu,
                      color='Aprox. KL',
                      fill='Aprox. KL'),alpha=0.8)+
    geom_point(aes(x=1:n,y=y_pred_mean_true,color='Verdadeiro',fill='Verdadeiro'),linetype='dashed',shape=4)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl_true,
                      ymax=y_pred_icu_true,
                      color='Verdadeiro',
                      fill='Verdadeiro'),alpha=0.8,linetype='dashed')+
    scale_x_continuous('Índice da amostra',limits=c(0,100))+
    scale_y_continuous('Valor observado')+
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    labs(title=TeX('Distribuição preditiva dos $y$'))+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}
vals_mean_true=data_mean%*%beta_true
vals_std_true=exp(-data_var%*%beta_true/2)

vals_mean_mean=data_mean%*%beta_mean
vals_var_mean=data_var%*%beta_mean

vals_mean_var=diag(data_mean%*%beta_var%*%t(data_mean))
vals_var_var=diag(data_var%*%beta_var%*%t(data_var))

vals_cov=diag(data_var%*%beta_var%*%t(data_mean))

mu0 <- vals_mean_mean + vals_cov
c0 <- exp(-vals_var_mean - vals_var_var / 2) / (vals_mean_var)
helper <- -3 + 3 * sqrt(1 + 2 * vals_var_var / 3)
# helper=Qt[2,2]
# print(c0)
# print(ft)
# print(Qt)
alpha <- 1 / helper
beta <- alpha * exp(-vals_var_mean - vals_var_var / 2)

nu <- 2 * alpha
sigma2 <- (beta / alpha) * (1 + 1 / c0)

y_pred_mean=(mu0-vals_mean_true)/vals_std_true
y_pred_icl=(qt(0.025,nu) * sqrt(sigma2) + mu0-vals_mean_true)/vals_std_true
y_pred_icu=(qt(0.975,nu) * sqrt(sigma2) + mu0-vals_mean_true)/vals_std_true

vals_mean_sample=data_mean%*%beta_sample
vals_var_sample=data_var%*%beta_sample

y_sample=(rnorm(prod(dim(vals_mean_sample))) %>% matrix(n,N))*exp(-vals_var_sample/2)+vals_mean_sample

y_pred_mean_true=(rowMeans(y_sample)-vals_mean_true)/vals_std_true
y_pred_icl_true=(apply(y_sample,1,function(x){quantile(x,0.025)})-vals_mean_true)/vals_std_true
y_pred_icu_true=(apply(y_sample,1,function(x){quantile(x,0.975)})-vals_mean_true)/vals_std_true


(ggplot()+
    geom_point(aes(x=1:n,y=(y-vals_mean_true)/vals_std_true,color='Observações',fill='Observações'))+
    geom_point(aes(x=1:n,y=y_pred_mean,color='Aprox. KL',fill='Aprox. KL'),linetype='dashed',shape=3)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl,
                      ymax=y_pred_icu,
                      color='Aprox. KL',
                      fill='Aprox. KL'),alpha=0.8)+
    geom_point(aes(x=1:n,y=y_pred_mean_true,color='Verdadeiro',fill='Verdadeiro'),linetype='dashed',shape=4)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl_true,
                      ymax=y_pred_icu_true,
                      color='Verdadeiro',
                      fill='Verdadeiro'),alpha=0.8,linetype='dashed')+
    scale_x_continuous('Índice da amostra',limits=c(0,100))+
    scale_y_continuous('Valor observado')+
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    labs(title=TeX('Distribuição preditiva dos $y$ padronizados'))+
    theme_bw()+
    theme(legend.position="bottom"))
```

Claramente há uma melhora substancial na qualidade da aproximação para um conjunto de dados mais extenso. Contudo, de forma análoga, temos uma redução substancial na qualidade do ajuste quanto tomamos um conjunto de dados pequeno. Para exemplificar esse comportamento, vamos exibir o mesmo ajuste feito anteriormente, mas agora para uma amostra de tamanho $20$.

```{r}
set.seed(13031998)

n=20

data_mean=data.frame(
  var0=1, # Intercepto média
  var1=rgamma(n,5,5/2),
  var2=rnorm(n,3,2),
  var3=0, # Intercepto variância
  var4=0
) %>% as.matrix

data_var=data.frame(
  var0=0, # Intercepto média
  var1=0,
  var2=0,
  var3=rep(1,n), # Intercepto variância
  var4=rpois(n,20)
) %>% as.matrix

beta_true=c(1,0.5,0.25,1,-0.04)
k=length(beta_true)

lin_pred_mean=as.matrix(data_mean)%*%beta_true
lin_pred_var=as.matrix(data_var)%*%beta_true

y=rnorm(n,lin_pred_mean,exp(-lin_pred_var/2))

# index=5
# plot(data_mean[,index]+data_var[,index],y)
```


```{r}
beta_mean=beta_mean_prior=rep(0,k)
beta_var=beta_var_prior=1*diag(k)

for(i in 1:n){
  FF=cbind(data_mean[i,],data_var[i,])
  x_mean=t(FF)%*%beta_mean
  Sigma=t(FF)%*%beta_var%*%FF
  mu1=x_mean[1]
  mu2=x_mean[2]
  s1=Sigma[1,1]
  s2=Sigma[2,2]
  rho=Sigma[1,2]/sqrt(s1*s2)
  val=rho*sqrt(s1/s2)

  x_var_adapt=(1-rho**2)*s1
  prec_prior=1/x_var_adapt
  f=function(x){
    x_mean_adapt=mu1+val*(log(x)-mu2)

    prec_obs=x
    x_var_update=1/(prec_obs+prec_prior)
    w1=1/(1+prec_prior/prec_obs)
    w2=1/(prec_obs/prec_prior+1)
    x_mean_update=w1*y[i]+w2*x_mean_adapt
    var_obs=1/prec_obs+x_var_adapt

    vals=dnorm(y[i],x_mean_adapt,sqrt(var_obs))*dlnorm(x,mu2,sqrt(s2))
    rbind(vals,
          log(x)*vals,
          (log(x)**2)*vals,
          x_mean_update*vals,
          (x_mean_update**2+x_var_update)*vals,
          (x_mean_update*log(x))*vals
    )

  }
  int_prev=cubintegrate(f,0,Inf,nVec=200,fDim=6)
  vals_prev=int_prev$integral
  c_prev=vals_prev[1]

  y_mean_prev=vals_prev[2]/c_prev
  y2_mean_prev=vals_prev[3]/c_prev
  y_s2_prev=y2_mean_prev-y_mean_prev**2

  x_mean_prev=vals_prev[4]/c_prev
  x2_mean_prev=vals_prev[5]/c_prev
  x_s2_prev=x2_mean_prev-x_mean_prev**2

  xy_mean_prev=vals_prev[6]/c_prev
  cov_prev=xy_mean_prev-x_mean_prev*y_mean_prev

  # sample=mvtnorm::rmvnorm(1000000,x_mean,Sigma)
  # mean((sample[,1]*sample[,2])*dnorm(y[i],sample[,1],exp(-sample[,2]/2)))/mean(dnorm(y[i],sample[,1],exp(-sample[,2]/2)))

  x_mean_new=c(x_mean_prev,y_mean_prev)
  Sigma_new=matrix(c(x_s2_prev,cov_prev,cov_prev,y_s2_prev),2,2)

  At <- beta_var %*% FF %*% ginv(Sigma)
  beta_mean <- beta_mean + At %*% (x_mean_new - x_mean)
  beta_var <- beta_var + At %*% (Sigma_new - Sigma) %*% t(At)
}
```


```{r}
log.like=function(beta){
  vals_mean=data_mean%*%beta
  vals_var=data_var%*%beta
  sum(dnorm(y,vals_mean,exp(-vals_var/2),log=TRUE))+dmvnorm(t(beta),beta_mean_prior,beta_var_prior,log=TRUE)
}

d_log_like=function(beta){calculus::derivative(function(beta){log.like(t(beta))},var=t(beta))}

mode=multiroot(d_log_like,start=beta_mean)$root

A=solve(-calculus::hessian(function(beta){log.like(t(beta))},var=t(mode)))
```

```{r}
beta_mean_ref=mode
beta_var_ref=A

N=5000
beta_sample=matrix(NA,k,N)
acep_list=matrix(0,k,N)
acep_p_list=matrix(0,k,N)
beta_unit=beta_mean_ref
cur_log_like=log.like(beta_unit)
beta_var_list=rep(NA,k)
for(j in 1:k){
  beta_var_list[j]=beta_var_ref[j,j]-t(beta_var_ref[j,-j])%*%ginv(beta_var_ref[-j,-j])%*%beta_var_ref[j,-j]
}
for(i in 1:N){
  for(j in 1:k){
    beta_mean_j=beta_mean_ref[j]+beta_var_ref[j,-j]%*%ginv(beta_var_ref[-j,-j])%*%(beta_unit[-j]-beta_mean_ref[-j])
    beta_var_j=beta_var_list[j]
    prop_j=rnorm(1,beta_mean_j,sqrt(beta_var_j))
    prop=beta_unit
    prop[j]=prop_j
    prop_log_like=log.like(prop)
    prop_log_val=dnorm(prop_j,beta_mean_j,sqrt(beta_var_j),log=TRUE)
    cur_log_val=dnorm(beta_unit[j],beta_mean_j,sqrt(beta_var_j),log=TRUE)
    acep=exp(prop_log_like-cur_log_like+cur_log_val-prop_log_val)
    acep_list[j,i]=acep
    if(runif(1)<acep){
      beta_unit[j]=prop_j
      cur_log_like=prop_log_like
      acep_list[j,i]=1
    }
  }
  beta_sample[,i]=beta_unit
}
beta_mean_true=rowMeans(beta_sample)
beta_var_true=var(t(beta_sample))

# mean(acep_list)

```



```{r}
data_plot=data.frame()
data_sample=pivot_longer(cbind(1:N,as.data.frame(t(beta_sample))),-1)
names(data_sample)=c('samp','index','value')
data_sample$index=substr(data_sample$index,2,2)
for(index in 1:k){
  x=seq(qnorm(1/N,beta_mean_true[index],sqrt(beta_var_true[index,index])),
        qnorm(1-1/N,beta_mean_true[index],sqrt(beta_var_true[index,index])),
        l=1000)
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,beta_mean[index],sqrt(beta_var[index,index])),
    label='Aprox. KL alt.',
    index=as.character(index)
  ))
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,mode[index],sqrt(A[index,index])),
    label='Aprox. Laplace',
    index=as.character(index)
  ))
  data_plot=rbind(data_plot,data.frame(
    x=x,
    fx=dnorm(x,beta_mean_true[index],sqrt(beta_var_true[index,index])),
    label='Aprox. KL orig.',
    index=as.character(index)
  ))

}

data_plot$index=as.factor(data_plot$index)
data_sample$index=as.factor(data_sample$index)
levels(data_sample$index)=levels(data_plot$index)=c('beta_0',
                                                    'beta_1',
                                                    'beta_2',
                                                    'beta_3',
                                                    'beta_4')
levels(data_sample$index)=levels(data_plot$index)=c('beta_0'=TeX('\\beta_0'),
                                                    'beta_1'=TeX('\\beta_1'),
                                                    'beta_2'=TeX('\\beta_2'),
                                                    'beta_3'=TeX('\\beta_3'),
                                                    'beta_4'=TeX('\\beta_4'))

(ggplot(data_plot %>% filter(label!='Aprox. Laplace'))+
    geom_histogram(aes(x=value,y=..density..),color='black',fill='#aaaaaa',bins=30,data=data_sample)+
    geom_line(aes(x=x,y=fx,color=label,linetype=label=='Aprox. KL orig.'))+
    scale_color_manual('',values=c('blue','red','green'))+
    scale_y_continuous(expand=c(0,0,0,0.5),limits=c(0,NA))+
    facet_wrap(.~index,scales='free')+
    guides(linetype='none')+
    theme_bw()+
    theme(legend.position="bottom"))
```

A princípio pode parecer que a perda de qualidade da aproximação não é tão grande com uma amostra de tamanho $20$, mas é importante que o leitor repare a escala do eixo $x$. Apesar das curvas parecerem próximas, a posteriori dos parâmetros está bem mais dispersa, de modo que as diferenças entre as distribuições correta e aproximada são de fato significativas. A redução da qualidade da aproximação é especialmente notável na distribuição preditiva de $y$:

  ```{r}

vals_mean_mean=data_mean%*%beta_mean
vals_var_mean=data_var%*%beta_mean

vals_mean_var=diag(data_mean%*%beta_var%*%t(data_mean))
vals_var_var=diag(data_var%*%beta_var%*%t(data_var))

vals_cov=diag(data_var%*%beta_var%*%t(data_mean))

mu0 <- vals_mean_mean + vals_cov
c0 <- exp(-vals_var_mean - vals_var_var / 2) / (vals_mean_var)
helper <- -3 + 3 * sqrt(1 + 2 * vals_var_var / 3)
# helper=Qt[2,2]
# print(c0)
# print(ft)
# print(Qt)
alpha <- 1 / helper
beta <- alpha * exp(-vals_var_mean - vals_var_var / 2)

nu <- 2 * alpha
sigma2 <- (beta / alpha) * (1 + 1 / c0)

y_pred_mean=mu0
y_pred_icl=qt(0.025,nu) * sqrt(sigma2) + mu0
y_pred_icu=qt(0.975,nu) * sqrt(sigma2) + mu0

vals_mean_sample=data_mean%*%beta_sample
vals_var_sample=data_var%*%beta_sample

y_sample=(rnorm(prod(dim(vals_mean_sample))) %>% matrix(n,N))*exp(-vals_var_sample/2)+vals_mean_sample

y_pred_mean_true=rowMeans(y_sample)
y_pred_icl_true=apply(y_sample,1,function(x){quantile(x,0.025)})
y_pred_icu_true=apply(y_sample,1,function(x){quantile(x,0.975)})


(ggplot()+
    geom_point(aes(x=1:n,y=y,color='Observações',fill='Observações'))+
    geom_point(aes(x=1:n,y=y_pred_mean,color='Aprox. KL',fill='Aprox. KL'),linetype='dashed',shape=3)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl,
                      ymax=y_pred_icu,
                      color='Aprox. KL',
                      fill='Aprox. KL'),alpha=0.8)+
    geom_point(aes(x=1:n,y=y_pred_mean_true,color='Verdadeiro',fill='Verdadeiro'),linetype='dashed',shape=4)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl_true,
                      ymax=y_pred_icu_true,
                      color='Verdadeiro',
                      fill='Verdadeiro'),alpha=0.8,linetype='dashed')+
    scale_x_continuous('Índice da amostra')+
    scale_y_continuous('Valor observado')+
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    labs(title=TeX('Distribuição preditiva dos $y$'))+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}
vals_mean_true=data_mean%*%beta_true
vals_std_true=exp(-data_var%*%beta_true/2)

vals_mean_mean=data_mean%*%beta_mean
vals_var_mean=data_var%*%beta_mean

vals_mean_var=diag(data_mean%*%beta_var%*%t(data_mean))
vals_var_var=diag(data_var%*%beta_var%*%t(data_var))

vals_cov=diag(data_var%*%beta_var%*%t(data_mean))

mu0 <- vals_mean_mean + vals_cov
c0 <- exp(-vals_var_mean - vals_var_var / 2) / (vals_mean_var)
helper <- -3 + 3 * sqrt(1 + 2 * vals_var_var / 3)
# helper=Qt[2,2]
# print(c0)
# print(ft)
# print(Qt)
alpha <- 1 / helper
beta <- alpha * exp(-vals_var_mean - vals_var_var / 2)

nu <- 2 * alpha
sigma2 <- (beta / alpha) * (1 + 1 / c0)

y_pred_mean=(mu0-vals_mean_true)/vals_std_true
y_pred_icl=(qt(0.025,nu) * sqrt(sigma2) + mu0-vals_mean_true)/vals_std_true
y_pred_icu=(qt(0.975,nu) * sqrt(sigma2) + mu0-vals_mean_true)/vals_std_true

vals_mean_sample=data_mean%*%beta_sample
vals_var_sample=data_var%*%beta_sample

y_sample=(rnorm(prod(dim(vals_mean_sample))) %>% matrix(n,N))*exp(-vals_var_sample/2)+vals_mean_sample

y_pred_mean_true=(rowMeans(y_sample)-vals_mean_true)/vals_std_true
y_pred_icl_true=(apply(y_sample,1,function(x){quantile(x,0.025)})-vals_mean_true)/vals_std_true
y_pred_icu_true=(apply(y_sample,1,function(x){quantile(x,0.975)})-vals_mean_true)/vals_std_true


(ggplot()+
    geom_point(aes(x=1:n,y=(y-vals_mean_true)/vals_std_true,color='Observações',fill='Observações'))+
    geom_point(aes(x=1:n,y=y_pred_mean,color='Aprox. KL',fill='Aprox. KL'),linetype='dashed',shape=3)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl,
                      ymax=y_pred_icu,
                      color='Aprox. KL',
                      fill='Aprox. KL'),alpha=0.8)+
    geom_point(aes(x=1:n,y=y_pred_mean_true,color='Verdadeiro',fill='Verdadeiro'),linetype='dashed',shape=4)+
    geom_errorbar(aes(x=1:n,
                      ymin=y_pred_icl_true,
                      ymax=y_pred_icu_true,
                      color='Verdadeiro',
                      fill='Verdadeiro'),alpha=0.8,linetype='dashed')+
    scale_x_continuous('Índice da amostra')+
    scale_y_continuous('Valor observado')+
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    labs(title=TeX('Distribuição preditiva dos $y$ padronizados'))+
    theme_bw()+
    theme(legend.position="bottom"))
```
