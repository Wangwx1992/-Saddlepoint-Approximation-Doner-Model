# 在 R 中执行以下命令选择清华镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("matlab")
#install.packages("polynom")
library(matlab)
library(polynom)
# function for solving 3-order equation
root3 <- function(a,b,c,d){
  xs <- polyroot(c(d,c,b,a))
  xx <- c()
  for (i in 1:3){
    if (abs(Im(xs[i]))<1.e-9){
      xx=c(xx,Re(xs[i]))
    }
  }
  x=xx[(xx>0&xx<1)]
  return(x)
} 

# function for finding mle estimates under null hypothesis
nullmle <- function(m){
  mp=apply(m,1,sum)
  mm=sum(mp)
  pihat = (mp[2]+2*mp[3])/mm/2
  rhohat = (4*mp[1]*mp[3]-mp[2]^2)/((mp[2]+2*mp[3])*(mp[2]+2*mp[1]))
  return (c(pihat,rhohat))
}

# loglikelihood
logl <- function(m,a,rho,d){
  aa=c(a,a+d)
  ll=sum(sum(m*log(rbind((1-aa)*(rho*aa-aa+1),2*aa*(1-aa)*(1-rho),rho*aa+(1-rho)*aa^2))))
}

#m <- matrix(c(15,3,13,14,9,21),nrow=3,ncol=2)
#nullmle(m)
#logl(m,0.4,0.9,0.2)

######wald 方法
# function for finding mle estimate of (pi1,pi2,rho)
findmle <- function(m){
  mc=apply(m,1,sum)
  mr=apply(m,2,sum)
  g=length(mr)
  rho=0
  nmle=nullmle(m)
  pihat=nmle[1]
  rho0 = nmle[2]
  pihat=ones(1,g)*pihat
  kk=0
  while (abs(rho-rho0)>1.e-6 && kk<200){
    kk=kk+1
    rho=c(rho0)
    for (ii in 1:g){
      a=(4*rho-2*rho^2-2)*mr[ii]; b=3*rho^2*mr[ii]-(c(5,6,7)*rho)%*%m[,ii]+c(2,3,4)%*%m[,ii]; c=(4*rho-rho^2)*mr[ii]-c(2*rho,1,2)%*%m[,ii]; d=-rho*c(0,1,1)%*%m[,ii];
      aaaa=root3(a,b,c,d)
      if (isempty(as.vector(aaaa))){
        rhohat=NaN
        break
      }
      if (length(as.vector(aaaa))==1){
        pihat[ii]=as.vector(aaaa)
      }
      else{
        prob=rbind((1-rho)*(1-as.vector(aaaa))^2+rho*(1-as.vector(aaaa)), 2*as.vector(aaaa)*(1-as.vector(aaaa))*(1-rho), rho*as.vector(aaaa)+(1-rho)*as.vector(aaaa)^2)
        aaaa=aaaa[sum(prob<0)==0]
      }
      if (length(as.vector(aaaa))==1){
        pihat[ii]=as.vector(aaaa);
      }
      else { 
        prob=rbind((1-rho)*(1-as.vector(aaaa))^2+rho*(1-as.vector(aaaa)), 2*as.vector(aaaa)*(1-as.vector(aaaa))*(1-rho), rho*as.vector(aaaa)+(1-rho)*as.vector(aaaa)^2)
        a=sum(log(prob)*kronecker(ones(1,dim(prob)[2]),m[,ii]))
        pihat[ii]=as.vector(aaaa)[a==max(a)]
      }
      ddrho= -1/(rho-1)^2*mc[2] - (pihat^2/(rho*pihat-pihat+1)^2)%*%m[1,]- ((pihat-1)^2/(rho+pihat-rho*pihat)^2)%*%m[3,]
      rho0=rho-1/ddrho *(1/(rho-1)*mc[2]+(pihat/(rho*pihat-pihat+1))%*%m[1,] - ((pihat-1)/(rho+pihat-rho*pihat))%*%m[3,])
    }
  }
  rhohat=rho
  return(c(pihat,rhohat))
}

#alpha =0.05
WaldCI <- function(m,alpha){
  m[m==0]=0.1
  est=findmle(m)
  diffest=est[2]- est[1];c=c(-1, 1, 0);a0=c(est[1],est[2]);rho=est[3];
  ddA = sum(m)*((-4*rho^2*a0^2 + 4*rho^2*a0-rho^2+6*rho*a0^2-6*rho*a0+2*rho+2*a0-2*a0^2)/(a0*(1-a0)*(rho+a0-rho*a0)*(rho*a0-a0+1)))
  ddArho =sum(m)*(rho*(2*a0-1)/((rho+a0-rho*a0)*(rho*a0-a0+1)))
  ddrho= sum(sum(m)*((a0*(1+rho)*(1-a0))/((1-rho)*(rho+a0-rho*a0)*(rho*a0-a0+1))))
  expInfoH0=rbind(cbind(diag(ddA),ddArho),cbind(t(ddArho),ddrho))
  ii=solve(expInfoH0)
  se=sqrt(t(c)%*%ii%*%c)
  lci=diffest-qnorm(1-alpha/2,0,1)*se
  uci=diffest+qnorm(1-alpha/2,0,1)*se
  CI=c(lci,uci)
  return(CI)
}


##似然比方法
###Fisher 矩阵
fishermatrix2 <- function(m,d,a,rho){
  ss=sum(m);m1=ss[1]; m2=ss[2];a0=c(a,a+d);
  #m01=m[1,1];m02=m[1,2];m11=m[2,1];m12=m[2,2];m21=m[3,1];m22=m[3,2];
  #a = est[1]
  #d = diffest
  #p01 = (1-a)*(rho*a-a+1); p11 = 2*a*(1-rho)*(1-a);p21 = a^2+rho*a*(1-a);
  #p02 = (1-(a+d))*(rho*(a+d)-(a+d)+1); p12 = 2*(a+d)*(1-rho)*(1-(a+d));p22 = (a+d)^2+rho*(a+d)*(1-(a+d));
  Idd = m2*((-4*rho^2*(a+d)^2 + 4*rho^2*(a+d)-rho^2+6*rho*(a+d)^2-6*rho*(a+d)+2*rho+2*(a+d)-2*(a+d)^2)/((a+d)*(1-(a+d))*(rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
  Iaa = m1*((-4*rho^2*(a)^2 + 4*rho^2*(a)-rho^2+6*rho*(a)^2-6*rho*(a)+2*rho+2*a-2*(a)^2)/((a)*(1-(a))*(rho+(a)-rho*(a))*(rho*(a)-(a)+1)))+m2*((-4*rho^2*(a+d)^2 + 4*rho^2*(a+d)-rho^2+6*rho*(a+d)^2-6*rho*(a+d)+2*rho+2*a0-2*(a+d)^2)/((a+d)*(1-(a+d))*(rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
  Irr = sum(sum(m)*((a0*(1+rho)*(1-a0))/((1-rho)*(rho+a0-rho*a0)*(rho*a0-a0+1))))
  Ida = Idd
  Idr = m2*(rho*(2*(a+d)-1)/((rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
  Iar = m1*(rho*(2*(a)-1)/((rho+(a)-rho*(a))*(rho*(a)-(a)+1))) + m2*(rho*(2*(a+d)-1)/((rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
  return(c(Idd,Iaa,Irr,Ida,Idr,Iar))
}

# find constrained mle
# 加载必要包
library(MASS) # 用于矩阵求逆

findMLE_constrained <- function(m,a0,rho0,d){
  rho=0;a=0;ss=sum(m); m1=ss[1]; m2=ss[2];
  m01=m[1,1];m02=m[1,2];m11=m[2,1];m12=m[2,2];m21=m[3,1];m22=m[3,2];
  kk=0;
  #d = d0
  while (abs(rho-rho0)>0.00001 && kk<200 ){
    kk=kk+1;
    rho=rho0;
    a=a0;
    b = a+d;
    #I = fishermatrix2(m,d,a,rho)
    #I11 = m1*((-4*rho^2*(a)^2 + 4*rho^2*(a)-rho^2+6*rho*(a)^2-6*rho*(a)+2*rho+2*a-2*(a)^2)/((a)*(1-(a))*(rho+(a)-rho*(a))*(rho*(a)-(a)+1)))+m2*((-4*rho^2*(a+d)^2 + 4*rho^2*(a+d)-rho^2+6*rho*(a+d)^2-6*rho*(a+d)+2*rho+2*a0-2*(a+d)^2)/((a+d)*(1-(a+d))*(rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
    #I12 = m1*(rho*(2*(a)-1)/((rho+(a)-rho*(a))*(rho*(a)-(a)+1))) + m2*(rho*(2*(a+d)-1)/((rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)))
    #I22 =  m1*((a*(1+rho)*(1-a))/((1-rho)*(rho+a-rho*a)*(rho*a-a+1))) + m2 *((b*(1+rho)*(1-b))/((1-rho)*(rho+b-rho*b)*(rho*b-b+1)))
    I11 = -m11*(2*a^2-2*a+1)/(a^2*(a-1)^2) - m21*(2*rho^2*a^2-2*rho^2*a+rho^2-4*rho*a^2+2*rho*a+2*a^2)/(a^2*(rho+a-rho*a)^2)- m01*(2*rho^2*a^2-2*rho^2*a+rho^2-4*rho*a^2+6*rho*a+2*a^2-2*rho+2)/((a-1)^2*(rho*a-a+1)^2) -
           m12*(2*b^2-2*b+1)/(b^2*(b-1)^2) - m22*(2*rho^2*b^2-2*rho^2*b+rho^2-4*rho*b^2+2*rho*b+2*b^2)/(b^2*(rho+b-rho*b)^2)- m02*(2*rho^2*b^2-2*rho^2*b+rho^2-4*rho*b^2+6*rho*b+2*b^2-2*rho+2)/((b-1)^2*(rho*b-b+1)^2);
    I12 = m01/(rho*a-a+1)^2-m21/(rho*a-a-rho)^2 + m02/(rho*b-b+1)^2-m22/(rho*b-b-rho)^2;
    I22 = -m11/(rho-1)^2-a^2*m01/(rho*a-a+1)^2-(a-1)^2*m21/(rho+a-rho*a)^2-m12/(rho-1)^2-b^2*m02/(rho*b-b+1)^2-(b-1)^2*m22/(rho+b-rho*b)^2;
    da = m11*(2*a-1)/(a*(a-1))+m21*(rho+2*a-2*rho*a)/(a*(rho+a-rho*a))-m01*(rho+2*a-2*rho*a-2)/((a-1)*(rho*a-a+1))+m12*(2*(a+d)-1)/((a+d)*((a+d)-1))+m22*(rho+2*(a+d)-2*rho*(a+d))/((a+d)*(rho+(a+d)-rho*(a+d)))-m01*(rho+2*(a+d)-2*rho*(a+d)-2)/(((a+d)-1)*(rho*(a+d)-(a+d)+1));
    drho= m11/(rho-1) - ((a-1)*m21)/(rho+a-rho*a) +a*m01/(rho*a-a+1)+ m12/(rho-1) - (((a+d)-1)*m22)/(rho+(a+d)-rho*(a+d)) +(a+d)*m01/(rho*(a+d)-(a+d)+1)
    arho = c(a,rho) - ginv(matrix(c(I11,I12,I12,I22),2,2)) %*% matrix(c(da,drho),2,1)
    a0 = arho[1,1]
    rho0 = arho[2,1]
  }
  pihat=a;
  rhohat=rho;
  return(c(pihat,rhohat))
}

# Profile likelihood CI
ProfileCI <- function(m,alpha){
  m[m==0]=0.1
  startest=findmle(m);a0=startest[1];R0=startest[3];d0=startest[2]-startest[1];
  d1=d0;flag=1;size=0.01;
  llnull=logl(m,a0,R0,d0);kk=0;
  while (size>1.e-5 & kk<10000){
    kk=kk+1
    d0=d0+flag*size;
    est=findMLE_constrained(m,a0,R0,d0);
    a0=est[1];
    R0=est[2];
    ll=logl(m,a0,R0,d0);
    if (2*flag*(llnull-ll)<flag*qchisq(1-alpha,1)){
      flag=flag;size=size;
    }
    else {
      flag=-flag;size=size*0.1;
    }
  }
  a0=startest[1];R0=startest[3];flag=-1;size=0.01; tt=0;
  while (size>1.e-5 & d1>0+size & tt<10000){
    tt=tt+1
    d1=d1+flag*size
    est=findMLE_constrained(m,a0,R0,d1)
    a0=est[1]
    R0=est[2]
    ll=logl(m,a0,R0,d1)
    if (2*(-flag)*(llnull-ll)<(-flag)*qchisq(1-alpha,1)){
      flag=flag;size=size;
    }
    else {
      flag=-flag;size=size*0.1;
    }
  }
  if (kk<10000 & tt<10000){
    LCI=d1;
    UCI=d0;
    CI=c(LCI,UCI);
  }
  else {
    CI=c(NA,NA)
  }
  return(CI)
}


# Score CI

# score function
SCORE <- function(m,a,rho,d){
  b=a+d
  ss=sum(m);m1=ss[1];m2=ss[2];
  m01=m[1,1];m02=m[1,2];m11=m[2,1];m12=m[2,2];m21=m[3,1];m22=m[3,2];
  I11 = m2*((-4*rho^2*(a+d)^2 + 4*rho^2*(a+d)-rho^2+6*rho*(a+d)^2-6*rho*(a+d)+2*rho+2*(a+d)-2*(a+d)^2)/((a+d)*(1-(a+d))*(rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1))) ;
  I22 = m1*((-4*rho^2*(a)^2 + 4*rho^2*(a)-rho^2+6*rho*(a)^2-6*rho*(a)+2*rho+2*a-2*(a)^2)/((a)*(1-(a))*(rho+(a)-rho*(a))*(rho*(a)-(a)+1)))+m2*((-4*rho^2*(a+d)^2 + 4*rho^2*(a+d)-rho^2+6*rho*(a+d)^2-6*rho*(a+d)+2*rho+2*(a+d)-2*(a+d)^2)/((a+d)*(1-(a+d))*(rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)));
  I33 = m1*((a*(1+rho)*(1-a))/((1-rho)*(rho+a-rho*a)*(rho*a-a+1))) + m2 *((b*(1+rho)*(1-b))/((1-rho)*(rho+b-rho*b)*(rho*b-b+1)));
  I12 = I11;
  I13 = m2*(rho*(2*(a+d)-1)/((rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)));
  I23 =  m1*(rho*(2*(a)-1)/((rho+(a)-rho*(a))*(rho*(a)-(a)+1))) + m2*(rho*(2*(a+d)-1)/((rho+(a+d)-rho*(a+d))*(rho*(a+d)-(a+d)+1)));
  invInfo = 1/(I11-matrix(c(I12,I13),1,2)%*%ginv(matrix(c(I22,I23,I23,I33),2,2))%*%matrix(c(I12,I13),2,1));
    #1/(I11-(I12^2*I33-I12*I13*I23+I13^2*I22)/(I22*I33-I23^2))
  dd = m12*(2*b-1)/(b*(b-1)) + m22*(rho+2*b-2*rho*b)/(b*(rho+b-rho*b)) - m02*(rho+2*b-2*rho*b-2)/((b-1)*(rho*b-b+1));
  chi=(dd^2)*invInfo;
}


ScoreCI <- function(m,alpha){
  m[m==0]=0.1
  startest=findmle(m);a0=startest[1];rho0=startest[3];d0=startest[2]- startest[1];
  d1=d0;flag=1;size=0.01;
  llnull=logl(m,a0,rho0,d0);kk=0;
  while (size>1.e-5 & kk<10000){
    kk=kk+1
    d0=d0+flag*size;
    est=findMLE_constrained(m,a0,rho0,d0);
    a0=est[1];
    rho0=est[2];
    chi=SCORE(m,a0,rho0,d0);
    if (flag*chi<flag*qchisq(1-alpha,1)){
      flag=flag;size=size;
    }
    else {
      flag=-flag;size=size*0.1;
    }
  }
  a0=startest[1];rho0=startest[3];flag=-1;size=0.01; tt=0;
  while (size>1.e-5 & d1>0+size & tt<10000){
    tt=tt+1
    d1=d1+flag*size
    est=findMLE_constrained(m,a0,rho0,d1)
    a0=est[1]
    rho0=est[2]
    chi=SCORE(m,a0,rho0,d1)
    if ((-flag)*chi<(-flag)*qchisq(1-alpha,1)){
      flag=flag;size=size;
    }
    else {
      flag=-flag;size=size*0.1;
    }
  }
  if (kk<10000 & tt<10000){
    LCI=d1;
    UCI=d0;
    CI=c(LCI,UCI);
  }
  else {
    CI=c(NA,NA)
  }
  return(CI)
}

# MoverCI using Wilson score Method

WilsonCI <- function(m,alpha){
  m[m==0]=0.1; ss=sum(m); z=qnorm(1-alpha/2);
  ntilda=2*ss+z^2; phat=(m[2,]+2*m[3,])/(2*ss); ptilda=(m[2,]+2*m[3,]+0.5*z^2)/ntilda
  plci=ptilda-z/ntilda*(sqrt(2*ss*phat*(1-phat)+z^2/4))
  puci=ptilda+z/ntilda*(sqrt(2*ss*phat*(1-phat)+z^2/4))
  est=findmle(m)
  rho=est[3]
  #pmle=est[1:2]
  ehat=1/(1+rho)
  LCI=phat[2]-phat[1] - sqrt((phat[2]-plci[2])^2/ehat+(puci[1]-phat[1])^2/ehat)
  UCI=phat[2]-phat[1] + sqrt((puci[2]-phat[2])^2/ehat+(phat[1]-plci[1])^2/ehat)
  CI=c(LCI,UCI)
  return(CI)
}
# MoverCI using Agresti-Coull Method
AgrestiCI <- function(m,alpha){
  m[m==0]=0.1; ss=sum(m); z=qnorm(1-alpha/2);
  ntilda=2*ss+4; phat=(m[2,]+2*m[3,])/(2*ss); ptilda=(m[2,]+2*m[3,]+2)/ntilda
  plci=ptilda-z*(sqrt(ptilda*(1-ptilda)/ntilda))
  puci=ptilda+z*(sqrt(ptilda*(1-ptilda)/ntilda))
  est=findmle(m)
  rho=est[3]
  #pmle=est[1:2]
  ehat=1/(1+rho)
  LCI=phat[2]-phat[1] - sqrt((phat[2]-plci[2])^2/ehat+(puci[1]-phat[1])^2/ehat)
  UCI=phat[2]-phat[1] + sqrt((puci[2]-phat[2])^2/ehat+(phat[1]-plci[1])^2/ehat)
  CI=c(LCI,UCI)
  return(CI)
}

#m <- matrix(c(15,3,13,14,9,21),nrow=3,ncol=2)
#ProfileCI(m,alpha=0.05)
#ScoreCI(m,alpha=0.05)
#WaldCI(m,alpha=0.05)
#WilsonCI(m,alpha=0.05)
#AgrestiCI(m,alpha=0.05)

# 计算群体概率
calculate_p <- function(pi_i, rho) {
  p0 <- (1-pi_i)*(rho*pi_i-pi_i+1)
  p1 <- 2*pi_i*(1-pi_i)*(1-rho)
  p2 <- rho*pi_i+(1-rho)*pi_i^2
  return(c(p0 = p0, p1 = p1, p2 = p2))
}

# 模拟生成数据
generate_data <- function(n_i, pi_i, rho) {
  p <- calculate_p(pi_i, rho)
  counts <- rmultinom(1, size = n_i, prob = p)
  return(c(n0 = counts[1], n1 = counts[2], n2 = counts[3]))
}
install.packages("foreach")
install.packages("doParallel")

library(foreach)
library(doParallel)

初始化并行环境 ---
n_cores <- detectCores() - 1  # 使用（总核心数-1）避免卡死系统
cl <- makeCluster(n_cores)
registerDoParallel(cl)

simulate_CI <- function(pi1, d, rho, n1, n2, N, alpha=0.05) {
  true_d <- d
  pi2 <- pi1 + d
  # 添加鞍点法结果收集
  methods <- c("Profile", "Score", "Wald", "Wilson", "Agresti")
  coverage <- setNames(numeric(length(methods)), methods)
  avg_length <- setNames(numeric(length(methods)), methods)
  counts <- setNames(numeric(length(methods)), methods)

  
  for (i in 1:N) {
    # 生成数据
    m1 <- generate_data(n1, pi1, rho)
    m2 <- generate_data(n2, pi2, rho)
    m <- cbind(m1, m2)
    
    # 计算各CI
    cis <- list(
      Profile = tryCatch(ProfileCI(m, alpha), error = function(e) c(NA, NA)),
      Score = tryCatch(ScoreCI(m, alpha), error = function(e) c(NA, NA)),
      Wald = tryCatch(WaldCI(m, alpha), error = function(e) c(NA, NA)),
      Wilson = tryCatch(WilsonCI(m, alpha), error = function(e) c(NA, NA)),
      Agresti = tryCatch(AgrestiCI(m, alpha), error = function(e) c(NA, NA))
    )
    
    # 更新统计量
    for (method in methods) {
      ci <- cis[[method]]
      if (any(is.na(ci))) next
      
      counts[[method]] <- counts[[method]] + 1
      lower <- ci[1]
      upper <- ci[2]
      coverage[[method]] <- coverage[[method]] + (lower <= true_d & upper >= true_d)
      avg_length[[method]] <- avg_length[[method]] + (upper - lower)
    }
    if (i %% 100 == 0) cat("已完成", i, "次模拟\n")
  }
  
  # 计算均值
  coverage <- coverage / counts
  avg_length <- avg_length / counts
  
  # 处理无效结果
  coverage[is.na(coverage)] <- 0
  avg_length[is.na(avg_length)] <- 0
  
  return(list(
    coverage = round(coverage, 3),
    avg_length = round(avg_length, 3),
    counts = counts
  ))
}

# 参数设置
params <- list(
  pi1 = 0.5,
  d = 0.2,
  rho = 0.9,
  n1 = 100,
  n2 = 100,
  N = 5000,
  alpha = 0.05
)

# 运行模拟
#set.seed(123)  # 确保结果可重复
results <- simulate_CI(
  params$pi1,
  params$d,
  params$rho,
  params$n1,
  params$n2,
  params$N,
  params$alpha
)

# 输出结果
cat("Simulation Results (N =", params$N, ")\n")
cat("True d:", params$d, "\n\n")
cat("Coverage Rates:\n")
print(results$coverage)
cat("\nAverage Lengths:\n")
print(results$avg_length)
cat("\nValid Counts:\n")
print(results$counts)

