
```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  fig.align='center'
)
Sys.setlocale("LC_TIME", "Portuguese")

library(tidyverse)
library(MASS)
library(rootSolve)
library(cubature)
library(mvtnorm)
library(latex2exp)
library(extraDistr)
```

## Caso Rayleight

```{r}
set.seed(13031998)

n=200

data=data.frame(
  var0=1, # Intercepto média
  var1=rgamma(n,5,5/2),
  var2=rbern(n,0.4)
) %>% as.matrix

beta_true=c(-1,-0.5,1)
k=length(beta_true)

lin_pred=as.matrix(data)%*%beta_true

y=rrayleigh(n,exp(lin_pred))

index=3
plot(y)
plot(data[,index],y)
hist(lin_pred)
```


```{r}
beta_mean=beta_mean_prior=rep(0,k)
beta_var=beta_var_prior=1*diag(k)

for(i in 1:n){
  FF=data[i,]
  x_mean=t(FF)%*%beta_mean
  Sigma=t(FF)%*%beta_var%*%FF

  f=function(x){
    vals=drayleigh(y[i],x)*dlnorm(x,x_mean,sqrt(Sigma))
    rbind(vals,
          log(x)*vals,
          (log(x)**2)*vals
    )
  }
  int_prev=cubintegrate(f,0,Inf,nVec=200,fDim=3)
  vals_prev=int_prev$integral
  c_prev=vals_prev[1]

  x_mean_prev=vals_prev[2]/c_prev
  x2_mean_prev=vals_prev[3]/c_prev
  x_s2_prev=y2_mean_prev-y_mean_prev**2

  # sample=mvtnorm::rmvnorm(1000000,x_mean,Sigma)
  # mean((sample[,1])*drayleigh(y[i],exp(sample[,1])))/mean(drayleigh(y[i],exp(sample[,1])))
  # mean((sample[,1]**2)*drayleigh(y[i],exp(sample[,1])))/mean(drayleigh(y[i],exp(sample[,1])))

  x_mean_new=x_mean_prev
  Sigma_new=x_s2_prev

  At <- beta_var %*% FF %*% ginv(Sigma)
  beta_mean <- beta_mean + At %*% (x_mean_new - x_mean)
  beta_var <- beta_var + At %*% (Sigma_new - Sigma) %*% t(At)
}
```


```{r}
log.like=function(beta){
  vals=data%*%beta
  sum(drayleigh(y,exp(vals),log=TRUE))+dmvnorm(t(beta),beta_mean_prior,beta_var_prior,log=TRUE)
}

d_log_like=function(beta){calculus::derivative(function(beta){log.like(t(beta))},var=t(beta))}

mode=multiroot(d_log_like,start=beta_true)$root

A=solve(-calculus::hessian(function(beta){log.like(t(beta))},var=t(mode)))
```

```{r}
beta_mean_ref=beta_mean
beta_var_ref=beta_var

N=500
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
mean(acep_list)
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
    scale_color_manual('',values=c('black','red','green'))+
    scale_y_continuous(expand=c(0,0,0,0.5),limits=c(0,NA))+
    facet_wrap(.~index,scales='free')+
    guides(linetype='none')+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}
beta_sample_alt=rmvnorm(N,beta_mean,beta_var)
vals_sample=data%*%t(beta_sample_alt)

y_sample=(rrayleigh(prod(dim(vals_sample))) %>% matrix(n,N))*exp(vals_sample)

y_pred_mean=rowMeans(y_sample)
y_pred_icl=apply(y_sample,1,function(x){quantile(x,0.025)})
y_pred_icu=apply(y_sample,1,function(x){quantile(x,0.975)})

vals_sample=data%*%beta_sample

y_sample=(rrayleigh(prod(dim(vals_sample))) %>% matrix(n,N))*exp(vals_sample)

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
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}
beta_sample_alt=rmvnorm(N,beta_mean,beta_var)
vals_sample=data%*%t(beta_sample_alt)

y_sample=(rrayleigh(prod(dim(vals_sample))) %>% matrix(n,N))*exp(vals_sample)

y_pred_mean=rowMeans(y_sample)/exp(lin_pred)
y_pred_icl=apply(y_sample,1,function(x){quantile(x,0.025)})/exp(lin_pred)
y_pred_icu=apply(y_sample,1,function(x){quantile(x,0.975)})/exp(lin_pred)

vals_sample=data%*%beta_sample

y_sample=(rrayleigh(prod(dim(vals_sample))) %>% matrix(n,N))*exp(vals_sample)

y_pred_mean_true=rowMeans(y_sample)/exp(lin_pred)
y_pred_icl_true=apply(y_sample,1,function(x){quantile(x,0.025)})/exp(lin_pred)
y_pred_icu_true=apply(y_sample,1,function(x){quantile(x,0.975)})/exp(lin_pred)


(ggplot()+
    geom_point(aes(x=1:n,y=y/exp(lin_pred),color='Observações',fill='Observações'))+
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
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    theme_bw()+
    theme(legend.position="bottom"))
```

## Caso Laplace Assimétrica

```{r}
set.seed(13031998)

n=100
sigmoid=function(x){1/(1+exp(-x))}

data_scale=data.frame(
  var0=1, # Intercepto média
  var1=rgamma(n,5,5/2),
  var2=0, # Intercepto variância
  var3=0
) %>% as.matrix

data_p=data.frame(
  var0=0, # Intercepto média
  var1=0,
  var2=rep(1,n), # Intercepto variância
  var3=rpois(n,40)
) %>% as.matrix

beta_true=c(2,-0.3,0.1*40,-0.1)
k=length(beta_true)

lin_pred_scale=as.matrix(data_scale)%*%beta_true
lin_pred_p=as.matrix(data_p)%*%beta_true

y=rep(NA,n)
for(i in 1:n){
  y[i]=ald::rALD(1,mu=0,sigma=exp(lin_pred_scale)[i,1],p=sigmoid(lin_pred_p)[i,1])
}
vec_min=function(x,y){
  pre_min=x-y
  (pre_min-abs(pre_min))/2+y
}
vec_max=function(x,y){
  -vec_min(-x,-y)
}

dALD=function(x,sigma=1,p=0.5,log=FALSE){
  const=p*(1-p)/sigma
  val1=-(p)*x/sigma
  val2=-(p-1)*x/sigma
  val=vec_min(val1,val2)
  if(log){
    return(log(const)+val)
  }else{
    return(const*exp(val))
  }
}

hist(exp(lin_pred_scale),breaks=30)
hist(sigmoid(lin_pred_p),breaks=20)
```


```{r}
beta_mean_prior=rep(0,k)
beta_var_prior=1*diag(k)

log.like=function(beta){
  vals_scale=data_scale%*%beta
  vals_p=data_p%*%beta
  sum(dALD(y,exp(vals_scale),sigmoid(vals_p),log=TRUE))+dmvnorm(t(beta),beta_mean_prior,beta_var_prior,log=TRUE)
}

d_log_like=function(beta){calculus::derivative(function(beta){log.like(t(beta))},var=t(beta))}

beta_mean=mode=multiroot(d_log_like,start=beta_true)$root
beta_var=A=solve(-calculus::hessian(function(beta){log.like(t(beta))},var=t(mode)))
```

```{r}
beta_mean_ref=beta_mean
beta_var_ref=beta_var

N=500
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
mean(acep_list)
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

(ggplot(data_plot %>% filter(label=='Aprox. KL orig.'))+
    geom_histogram(aes(x=value,y=..density..),color='black',fill='#aaaaaa',bins=30,data=data_sample)+
    geom_line(aes(x=x,y=fx,color=label,linetype=label=='Aprox. KL orig.'))+
    scale_color_manual('',values=c('black','red','green'))+
    scale_y_continuous(expand=c(0,0,0,0.5),limits=c(0,NA))+
    facet_wrap(.~index,scales='free')+
    guides(linetype='none')+
    theme_bw()+
    theme(legend.position="bottom"))
```

```{r}
beta_sample_alt=rmvnorm(N,beta_mean,beta_var)
vals_scale_sample=data_scale%*%t(beta_sample_alt)
vals_p_sample=data_p%*%t(beta_sample_alt)

y_sample=matrix(NA,n,N)
for(i in 1:n){
  for(j in 1:N){
    y_sample[i,j]=ald::rALD(1,0,exp(vals_scale_sample[i,j]),sigmoid(vals_p_sample[i,j]))
  }
}

y_pred_mean=rowMeans(y_sample)
y_pred_icl=apply(y_sample,1,function(x){quantile(x,0.025)})
y_pred_icu=apply(y_sample,1,function(x){quantile(x,0.975)})


vals_scale_sample=data_scale%*%beta_sample
vals_p_sample=data_p%*%beta_sample

y_sample=matrix(NA,n,N)
for(i in 1:n){
  for(j in 1:N){
    y_sample[i,j]=ald::rALD(1,0,exp(vals_scale_sample[i,j]),sigmoid(vals_p_sample[i,j]))
  }
}

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
    scale_color_manual('',values=c('#4444ff','black','#ff4444'))+
    scale_fill_manual('',values=c('#4444ff','black','#ff4444'))+
    theme_bw()+
    theme(legend.position="bottom"))
```
