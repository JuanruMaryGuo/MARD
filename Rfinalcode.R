library(distributions3)
library(stats)
library(resample)
library(meta)
library(parallel)
library(MASS)
library(xtable)
library(prevalence)


quadratic_roots <- function(aa = 1, bb = 2, cc = 1)
{
  # calculate the discriminant
  disc <- bb^2 - 4 * aa * cc
  
  if (disc == 0) # if1, single root
  {
    return(-bb / (2 * aa))
  }   
  else
  {
    if (disc > 0) # if2, two roots
    {
      root1 <- (-bb - sqrt(disc)) / (2 * aa)
      root2 <- (-bb + sqrt(disc)) / (2 * aa)
      return(c(root1, root2))
    }
    else # no real roots, return NA
    {
      return(NA)
    }   
  }   
}
 

find_mu_sigma<-function(mean,lci){
  m = qnorm(mean)
  l = qnorm(lci)
  answer = quadratic_roots(m*m-qnorm(0.975)**2,-2*qnorm(0.975)*l,m*m-l*l)
  sigma = answer[answer>0]
  mu = qnorm(lci)+qnorm(0.975)*sigma
  return(cbind(mu,sigma))
}



for (dz in c(1,2)){
  for (rho in c(0.3,0.6)){
    set.seed(620)
    mean_e = mean_es[dz]
    lower_e = lower_es[dz]
    mean_c = mean_cs[dz]
    lower_c = lower_cs[dz]
    var_e = tryCatch(find_mu_sigma(mean_e,lower_e)[2,2],error=function(e){find_mu_sigma(mean_e,lower_e)[1,2]})**2
    var_c = tryCatch(find_mu_sigma(mean_c,lower_c)[2,2],error=function(e){find_mu_sigma(mean_c,lower_c)[1,2]})**2
    Var <- matrix(c(var_e,sqrt(var_e*var_c)*rho,sqrt(var_e*var_c)*rho,var_c),2,2)
    mean = c(tryCatch(find_mu_sigma(mean_e,lower_e)[2,1],error=function(e){find_mu_sigma(mean_e,lower_e)[1,1]})
             ,tryCatch(find_mu_sigma(mean_c,lower_c)[2,1],error=function(e){find_mu_sigma(mean_c,lower_c)[1,1]}))
    list=pnorm(mvrnorm(Num, mean, Var))
    print(format(round(mean, 10), nsmall = 10))
    print(format(round(sqrt(Var), 10), nsmall = 10))
    print(colMeans(list))
  }
}



# calculate risk difference without double zero events
meta_withoutDZE_function<-function(databin){
  a = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "DL", incr = 0, data = databin)
  a1 = tryCatch(a[["TE.fixed"]], error=function(e){NA})
  a2 = tryCatch(a[["lower.fixed"]], error=function(e){NA})
  a3 = tryCatch(a[["upper.fixed"]], error=function(e){NA})
  b = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "MH", method.tau = "DL", incr = 0, data = databin)
  b1 = tryCatch(b[["TE.fixed"]], error=function(e){NA})
  b2 = tryCatch(b[["lower.fixed"]], error=function(e){NA})
  b3 = tryCatch(b[["upper.fixed"]], error=function(e){NA})
  c = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "DL", incr = 0, data = databin)
  c1 = tryCatch(c[["TE.random"]], error=function(e){NA})
  c2 = tryCatch(c[["lower.random"]], error=function(e){NA})
  c3 = tryCatch(c[["upper.random"]], error=function(e){NA})
  d = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HE", incr = 0, data = databin)
  d1 = tryCatch(d[["TE.random"]], error=function(e){NA})
  d2 = tryCatch(d[["lower.random"]], error=function(e){NA})
  d3 = tryCatch(d[["upper.random"]], error=function(e){NA})
  e = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "PM", incr = 0, data = databin)
  e1 = tryCatch(e[["TE.random"]], error=function(e){NA})
  e2 = tryCatch(e[["lower.random"]], error=function(e){NA})
  e3 = tryCatch(e[["upper.random"]], error=function(e){NA})
  f = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "ML", incr = 0, data = databin)
  f1 = tryCatch(f[["TE.random"]], error=function(e){NA})
  f2 = tryCatch(f[["lower.random"]], error=function(e){NA})
  f3 = tryCatch(f[["upper.random"]], error=function(e){NA})
  g = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "REML", incr = 0, data = databin)
  g1 = tryCatch(g[["TE.random"]], error=function(e){NA})
  g2 = tryCatch(g[["lower.random"]], error=function(e){NA})
  g3 = tryCatch(g[["upper.random"]], error=function(e){NA})
  h = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HS", incr = 0, data = databin)
  h1 = tryCatch(h[["TE.random"]], error=function(e){NA})
  h2 = tryCatch(h[["lower.random"]], error=function(e){NA})
  h3 = tryCatch(h[["upper.random"]], error=function(e){NA})
  i = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "SJ", incr = 0, data = databin)
  i1 = tryCatch(i[["TE.random"]], error=function(e){NA})
  i2 = tryCatch(i[["lower.random"]], error=function(e){NA})
  i3 = tryCatch(i[["upper.random"]], error=function(e){NA})
  j = metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "EB", incr = 0, data = databin)
  j1 = tryCatch(j[["TE.random"]], error=function(e){NA})
  j2 = tryCatch(j[["lower.random"]], error=function(e){NA})
  j3 = tryCatch(j[["upper.random"]], error=function(e){NA})
  return(c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,
           a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,
           a3,b3,c3,d3,e3,f3,g3,h3,i3,j3)
  )
}

# calculate risk diffference with double zero events
meta_withDZE_function<-function(databin,databin_deletezeros){
  a = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", incr = 0, data = databin_deletezeros), error=function(e){NA})
  a1 = tryCatch(a[["TE.fixed"]], error=function(e){NA})
  a2 = tryCatch(a[["lower.fixed"]], error=function(e){NA})
  a3 = tryCatch(a[["upper.fixed"]], error=function(e){NA})
  b = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", incr = 0.5,  data = databin), error=function(e){NA})
  b1 = tryCatch(b[["TE.fixed"]], error=function(e){NA})
  b2 = tryCatch(b[["lower.fixed"]], error=function(e){NA})
  b3 = tryCatch(b[["upper.fixed"]], error=function(e){NA})
  c = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", incr = "TACC", data = databin), error=function(e){NA})
  c1 = tryCatch(c[["TE.fixed"]], error=function(e){NA})
  c2 = tryCatch(c[["lower.fixed"]], error=function(e){NA})
  c3 = tryCatch(c[["upper.fixed"]], error=function(e){NA})
  d = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "MH", data = databin), error=function(e){NA})
  d1 = tryCatch(d[["TE.fixed"]], error=function(e){NA})
  d2 = tryCatch(d[["lower.fixed"]], error=function(e){NA})
  d3 = tryCatch(d[["upper.fixed"]], error=function(e){NA})
  e = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "DL", incr = 0, data = databin_deletezeros), error=function(e){NA})
  e1 = tryCatch(e[["TE.random"]], error=function(e){NA})
  e2 = tryCatch(e[["lower.random"]], error=function(e){NA})
  e3 = tryCatch(e[["upper.random"]], error=function(e){NA})
  f = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HE", incr = 0, data = databin_deletezeros), error=function(e){NA})
  f1 = tryCatch(f[["TE.random"]], error=function(e){NA})
  f2 = tryCatch(f[["lower.random"]], error=function(e){NA})
  f3 = tryCatch(f[["upper.random"]], error=function(e){NA})
  g = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "PM", incr = 0, data = databin_deletezeros), error=function(e){NA})
  g1 = tryCatch(g[["TE.random"]], error=function(e){NA})
  g2 = tryCatch(g[["lower.random"]], error=function(e){NA})
  g3 = tryCatch(g[["upper.random"]], error=function(e){NA})
  h = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "ML", incr = 0, data = databin_deletezeros), error=function(e){NA})
  h1 = tryCatch(h[["TE.random"]], error=function(e){NA})
  h2 = tryCatch(h[["lower.random"]], error=function(e){NA})
  h3 = tryCatch(h[["upper.random"]], error=function(e){NA})
  i = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "REML", incr = 0, data = databin_deletezeros), error=function(e){NA})
  i1 = tryCatch(i[["TE.random"]], error=function(e){NA})
  i2 = tryCatch(i[["lower.random"]], error=function(e){NA})
  i3 = tryCatch(i[["upper.random"]], error=function(e){NA})
  j = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HS", incr = 0, data = databin_deletezeros), error=function(e){NA})
  j1 = tryCatch(j[["TE.random"]], error=function(e){NA})
  j2 = tryCatch(j[["lower.random"]], error=function(e){NA})
  j3 = tryCatch(j[["upper.random"]], error=function(e){NA})
  k = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "SJ", incr = 0, data = databin_deletezeros), error=function(e){NA})
  k1 = tryCatch(k[["TE.random"]], error=function(e){NA})
  k2 = tryCatch(k[["lower.random"]], error=function(e){NA})
  k3 = tryCatch(k[["upper.random"]], error=function(e){NA})
  l = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "EB", incr = 0, data = databin_deletezeros), error=function(e){NA})
  l1 = tryCatch(l[["TE.random"]], error=function(e){NA})
  l2 = tryCatch(l[["lower.random"]], error=function(e){NA})
  l3 = tryCatch(l[["upper.random"]], error=function(e){NA})
  m = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "DL", incr = 0.5, data = databin), error=function(e){NA})
  m1 = tryCatch(m[["TE.random"]], error=function(e){NA})
  m2 = tryCatch(m[["lower.random"]], error=function(e){NA})
  m3 = tryCatch(m[["upper.random"]], error=function(e){NA})
  n = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HE", incr = 0.5, data = databin), error=function(e){NA})
  n1 = tryCatch(n[["TE.random"]], error=function(e){NA})
  n2 = tryCatch(n[["lower.random"]], error=function(e){NA})
  n3 = tryCatch(n[["upper.random"]], error=function(e){NA})
  o = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "PM", incr = 0.5, data = databin), error=function(e){NA})
  o1 = tryCatch(o[["TE.random"]], error=function(e){NA})
  o2 = tryCatch(o[["lower.random"]], error=function(e){NA})
  o3 = tryCatch(o[["upper.random"]], error=function(e){NA})
  p = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "ML", incr = 0.5, data = databin), error=function(e){NA})
  p1 = tryCatch(p[["TE.random"]], error=function(e){NA})
  p2 = tryCatch(p[["lower.random"]], error=function(e){NA})
  p3 = tryCatch(p[["upper.random"]], error=function(e){NA})
  q = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "REML", incr = 0.5, data = databin), error=function(e){NA})
  q1 = tryCatch(q[["TE.random"]], error=function(e){NA})
  q2 = tryCatch(q[["lower.random"]], error=function(e){NA})
  q3 = tryCatch(q[["upper.random"]], error=function(e){NA})
  r = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HS", incr = 0.5, data = databin), error=function(e){NA})
  r1 = tryCatch(r[["TE.random"]], error=function(e){NA})
  r2 = tryCatch(r[["lower.random"]], error=function(e){NA})
  r3 = tryCatch(r[["upper.random"]], error=function(e){NA})
  s = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "SJ", incr = 0.5, data = databin), error=function(e){NA})
  s1 = tryCatch(s[["TE.random"]], error=function(e){NA})
  s2 = tryCatch(s[["lower.random"]], error=function(e){NA})
  s3 = tryCatch(s[["upper.random"]], error=function(e){NA})
  t = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "EB", incr = 0.5, data = databin), error=function(e){NA})
  t1 = tryCatch(t[["TE.random"]], error=function(e){NA})
  t2 = tryCatch(t[["lower.random"]], error=function(e){NA})
  t3 = tryCatch(t[["upper.random"]], error=function(e){NA})
  u = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "DL", incr = "TACC", data = databin), error=function(e){NA})
  u1 = tryCatch(u[["TE.random"]], error=function(e){NA})
  u2 = tryCatch(u[["lower.random"]], error=function(e){NA})
  u3 = tryCatch(u[["upper.random"]], error=function(e){NA})
  v = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HE", incr = "TACC", data = databin), error=function(e){NA})
  v1 = tryCatch(v[["TE.random"]], error=function(e){NA})
  v2 = tryCatch(v[["lower.random"]], error=function(e){NA})
  v3 = tryCatch(v[["upper.random"]], error=function(e){NA})
  w = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "PM", incr = "TACC", data = databin), error=function(e){NA})
  w1 = tryCatch(w[["TE.random"]], error=function(e){NA})
  w2 = tryCatch(w[["lower.random"]], error=function(e){NA})
  w3 = tryCatch(w[["upper.random"]], error=function(e){NA})
  x = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "ML", incr = "TACC", data = databin), error=function(e){NA})
  x1 = tryCatch(x[["TE.random"]], error=function(e){NA})
  x2 = tryCatch(x[["lower.random"]], error=function(e){NA})
  x3 = tryCatch(x[["upper.random"]], error=function(e){NA})
  y = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "REML", incr = "TACC", data = databin), error=function(e){NA})
  y1 = tryCatch(y[["TE.random"]], error=function(e){NA})
  y2 = tryCatch(y[["lower.random"]], error=function(e){NA})
  y3 = tryCatch(y[["upper.random"]], error=function(e){NA})
  z = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "HS", incr = "TACC", data = databin), error=function(e){NA})
  z1 = tryCatch(z[["TE.random"]], error=function(e){NA})
  z2 = tryCatch(z[["lower.random"]], error=function(e){NA})
  z3 = tryCatch(z[["upper.random"]], error=function(e){NA})
  aa = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "SJ", incr = "TACC", data = databin), error=function(e){NA})
  aa1 = tryCatch(aa[["TE.random"]], error=function(e){NA})
  aa2 = tryCatch(aa[["lower.random"]], error=function(e){NA})
  aa3 = tryCatch(aa[["upper.random"]], error=function(e){NA})
  ab = tryCatch(metabin(Ee, Ne, Ec, Nc, sm = "RD", method = "I", method.tau = "EB", incr = "TACC", data = databin), error=function(e){NA})
  ab1 = tryCatch(ab[["TE.random"]], error=function(e){NA})
  ab2 = tryCatch(ab[["lower.random"]], error=function(e){NA})
  ab3 = tryCatch(ab[["upper.random"]], error=function(e){NA})
  return(c(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1,q1,r1,s1,t1,u1,v1,q1,x1,y1,z1,aa1,ab1,
           a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,p2,q2,r2,s2,t2,u2,v2,q2,x2,y2,z2,aa2,ab2,
           a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,k3,l3,m3,n3,o3,p3,q3,r3,s3,t3,u3,v3,q3,x3,y3,z3,aa3,ab3)
  )
}


# calculate alpha & beta from mean and variance
alpha_beta_function<-function(m,v){
  alpha = m*(((m*(1-m))/v)-1)
  beta = ((m*(1-m))/v-1)*(1-m)
  return(c(alpha,beta))
}

# sample bivariate-beta
# sample beta
sample_beta<-function(n,alpha,beta){
  X <- Beta(alpha, beta)
  return(random(X, n))
}

# calculate g
g_part_function<-function(p1,p2,w,alpha1,alpha2,beta1,beta2){
  mu1 = alpha1/(alpha1+beta1)
  mu2 = alpha2/(alpha2+beta2)
  return((1+w*(p1-mu1)*(p2-mu2)))
}

# caluclate w
w_function<-function(rho,alpha1,alpha2,beta1,beta2){
  var1 = alpha1*beta1/((alpha1+beta1)^2*(alpha1+beta1+1))
  var2 = alpha2*beta2/((alpha2+beta2)^2*(alpha2+beta2+1))
  return(rho/sqrt(var1*var2))
}

# calculate M
M_function<-function(w){
  return(1+w)
}

# sample bivariate-beta
sample_bibeta<-function(n,alpha1,alpha2,beta1,beta2,rho){
  w = w_function(rho,alpha1,alpha2,beta1,beta2)
  M = M_function(w)
  num = 0
  while (num<n) {
    p1 = sample_beta(n,alpha1,beta1)
    p2 = sample_beta(n,alpha2,beta2)
    g_part = g_part_function(p1,p2,w,alpha1,alpha2,beta1,beta2)
    flag = g_part/M
    u = runif(n, min = 0, max = 1)
    keep = (u<=flag)
    p1 = p1 * keep
    p2 = p2 * keep
    p1 = p1[p1!=0]
    p2 = p2[p2!=0]
    if (num == 0){
      samplingp1 = p1
      samplingp2 = p2
    } else {
      samplingp1 = append(samplingp1,p1)
      samplingp2 = append(samplingp2,p2)
    }
    num = num+sum(keep)
  }
  samplingp1 = head(samplingp1,n)
  samplingp2 = head(samplingp2,n)
  x_name <- "p1"
  y_name <- "p2"
  output <- data.frame(samplingp1,samplingp2)
  names(output) <- c(x_name,y_name)
  return(output)
}


Num = 5000# Sampling Size
ncol = c(10,28,28) * 3 # without double zero events : 10 methods; with double zero events : 28 methods 
mean_es = c(0.2,0.01,0.08) # means of experimental group
lower_es =  c(0.1,0.005,0.0001) # variance of experimental group
mean_cs = c(0.4,0.03,0.1) # means of control group
lower_cs =  c(0.2,0.01,0.0005) # variance of control group
N=100


for (dz in c(2,3)){
  for (rho in c(0.3,0.6)){
    set.seed(0)
    mean_e = mean_es[dz]
    lower_e = lower_es[dz]
    mean_c = mean_cs[dz]
    lower_c = lower_cs[dz]
    var_e = tryCatch(find_mu_sigma(mean_e,lower_e)[2,2],error=function(e){find_mu_sigma(mean_e,lower_e)[1,2]})**2
    var_c = tryCatch(find_mu_sigma(mean_c,lower_c)[2,2],error=function(e){find_mu_sigma(mean_c,lower_c)[1,2]})**2
    Var <- matrix(c(var_e,sqrt(var_e*var_c)*rho,sqrt(var_e*var_c)*rho,var_c),2,2)
    mean = c(tryCatch(find_mu_sigma(mean_e,lower_e)[2,1],error=function(e){find_mu_sigma(mean_e,lower_e)[1,1]})
             ,tryCatch(find_mu_sigma(mean_c,lower_c)[2,1],error=function(e){find_mu_sigma(mean_c,lower_c)[1,1]}))
    result <- mclapply(c(10,30),function(size){
      list=pnorm(mvrnorm(Num*size, mean,Var))
      Pe = matrix(list[,1], ncol = size)
      Pc = matrix(list[,2], ncol = size)
      Ne = rep(N,size)
      Nc = rep(N,size)
      result <- as.data.frame(matrix(0, ncol = ncol[dz], nrow = Num))
      for (study in c(1:Num)){
        Ee = rbinom(size, N, Pe[(study-1)*size+1:study*size])
        Ec = rbinom(size, N, Pc[(study-1)*size+1:study*size])
        data = data.frame("Ee" = Ee, "Ne" = Ne, "Ec" = Ec, "Nc" = Nc)
        if (dz==1){
          result[study,] = meta_withoutDZE_function(data)
        }else{
          #print(data)
          row_sub = apply(data, 1, function(row) all(row !=0 ))
          data_deletezeros = data[row_sub,]
          #if(nrow(data_deletezeros)!=0){
          #  result[study,] = meta_withDZE_function(data,data_deletezeros)
          #}else{
          result[study,] = meta_withDZE_function(data,data_deletezeros)
          #t = t + 1
          #}
        }
      }
      write.csv(result,paste("C:\\Users\\mary\\Desktop\\metacodes\\N",dz,rho,size,".csv",sep = ""))
    }
    ,mc.cores=1) # depend on how many cores your pc has
  }
}



for (dz in c(1)){
  for (rho in c(0.3,0.6)){
    set.seed(1)
    mean_e = mean_es[dz]
    lower_e = lower_es[dz]
    mean_c = mean_cs[dz]
    lower_c = lower_cs[dz]
    var_e = tryCatch(find_mu_sigma(mean_e,lower_e)[2,2],error=function(e){find_mu_sigma(mean_e,lower_e)[1,2]})**2
    var_c = tryCatch(find_mu_sigma(mean_c,lower_c)[2,2],error=function(e){find_mu_sigma(mean_c,lower_c)[1,2]})**2
    Var <- matrix(c(var_e,sqrt(var_e*var_c)*rho,sqrt(var_e*var_c)*rho,var_c),2,2)
    mean = c(tryCatch(find_mu_sigma(mean_e,lower_e)[2,1],error=function(e){find_mu_sigma(mean_e,lower_e)[1,1]})
             ,tryCatch(find_mu_sigma(mean_c,lower_c)[2,1],error=function(e){find_mu_sigma(mean_c,lower_c)[1,1]}))
    result <- mclapply(c(10,30),function(size){
      list=pnorm(mvrnorm(Num*size, mean,Var))
      Pe = matrix(list[,1], ncol = size)
      Pc = matrix(list[,2], ncol = size)
      Ne = rep(N,size)
      Nc = rep(N,size)
      result <- as.data.frame(matrix(0, ncol = ncol[dz], nrow = Num))
      for (study in c(1:Num)){
        Ee = rbinom(size, N, Pe[(study-1)*size+1:study*size])
        Ec = rbinom(size, N, Pc[(study-1)*size+1:study*size])
        data = data.frame("Ee" = Ee, "Ne" = Ne, "Ec" = Ec, "Nc" = Nc)
        if (dz==1){
          result[study,] = meta_withoutDZE_function(data)
        }else{
          #print(data)
          row_sub = apply(data, 1, function(row) all(row !=0 ))
          data_deletezeros = data[row_sub,]
          #if(nrow(data_deletezeros)!=0){
          #  result[study,] = meta_withDZE_function(data,data_deletezeros)
          #}else{
          result[study,] = meta_withDZE_function(data,data_deletezeros)
          #t = t + 1
          #}
        }
      }
      write.csv(result,paste("C:\\Users\\mary\\Desktop\\metacodes\\N1",dz,rho,size,".csv",sep = ""))
    }
    ,mc.cores=1) # depend on how many cores your pc has
  }
}


for (dz in c(1,2,3)){
  for (rho in c(0.3,0.6)){
    set.seed(0)
    mean_e = mean_es[dz]
    lower_e = lower_es[dz]
    mean_c = mean_cs[dz]
    lower_c = lower_cs[dz]
    alpha_e = betaExpert(mean_e,lower_e, p = 0.95, method = "mean")[["alpha"]]
    beta_e = betaExpert(mean_e,lower_e, p = 0.95, method = "mean")[["beta"]]
    alpha_c = betaExpert(mean_c,lower_c, p = 0.95, method = "mean")[["alpha"]]
    beta_c = betaExpert(mean_c,lower_c, p = 0.95, method = "mean")[["beta"]]
    print(alpha_e)
    print(alpha_c)
    print(beta_e)
    print(beta_c)
  }
}


for (dz in c(1,2,3)){
  for (rho in c(0.3,0.6)){
    set.seed(0)
    mean_e = mean_es[dz]
    lower_e = lower_es[dz]
    mean_c = mean_cs[dz]
    lower_c = lower_cs[dz]
    alpha_e = betaExpert(mean_e,lower_e, p = 0.95, method = "mean")[["alpha"]]
    beta_e = betaExpert(mean_e,lower_e, p = 0.95, method = "mean")[["beta"]]
    alpha_c = betaExpert(mean_c,lower_c, p = 0.95, method = "mean")[["alpha"]]
    beta_c = betaExpert(mean_c,lower_c, p = 0.95, method = "mean")[["beta"]]
    result <- mclapply(c(10,30),function(size){
      list=sample_bibeta(Num*size,alpha_e,alpha_c,beta_e,beta_c,rho)
      Pe = matrix(list$p1, ncol = size)
      Pc = matrix(list$p2, ncol = size)
      Ne = rep(N,size)
      Nc = rep(N,size)
      result <- as.data.frame(matrix(0, ncol = ncol[dz], nrow = Num))
      for (study in c(1:Num)){
        Ee = rbinom(size, N, Pe[(study-1)*size+1:study*size])
        Ec = rbinom(size, N, Pc[(study-1)*size+1:study*size])
        data = data.frame("Ee" = Ee, "Ne" = Ne, "Ec" = Ec, "Nc" = Nc)
        if (dz==1){
          result[study,] = meta_withoutDZE_function(data)
        }else{
          row_sub = apply(data, 1, function(row) all(row !=0 ))
          data_deletezeros = data[row_sub,]
          result[study,] = meta_withDZE_function(data,data_deletezeros)
        }
      }
      result1<<-result
      write.csv(result,paste("C:\\Users\\mary\\Desktop\\metacodes\\R",dz,rho,size,".csv",sep = ""))
    }
    ,mc.cores=1) # depend on how many cores your pc has
  }
}

