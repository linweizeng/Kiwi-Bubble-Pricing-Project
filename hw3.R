library('data.table')
library('reshape')

#2[A]A logit demand system
#data preparation
kiwi<-read.csv('kiwi_bubbles.csv',stringsAsFactors = FALSE)
kiwi$choice<-rep(0,nrow(kiwi))
for (i in 1:nrow(kiwi)) {
  if(kiwi[i,'purchase_KR']==1){kiwi[i,'choice']<-1}
  if(kiwi[i,'purchase_KB']==1){kiwi[i,'choice']<-2}
  if(kiwi[i,'purchase_MB']==1){kiwi[i,'choice']<-3}
  kiwi$price[i]<-unique(paste(kiwi$price_KR[i],kiwi$price_KB[i],kiwi$price_MB[i]))  
}
price<-unique(paste(kiwi$price_KR,kiwi$price_KB,kiwi$price_MB))
price1<-strsplit(price,' ')
product<-data.frame(scenario=paste0('s',1:68),bubble.1=0,kiwi.1=1,price.1=0,bubble.2=1,kiwi.2=1,price.2=0,bubble.3=1,kiwi.3=0,price.3=0)
for (i in 1:nrow(product)) {
  product[i,'price.1']<-price1[[i]][1]
  product[i,'price.2']<-price1[[i]][2]
  product[i,'price.3']<-price1[[i]][3]
  product$price[[i]]<-price[[i]]
}
kiwi$scenario<-product$scenario[match(kiwi$price,product$price)]
kiwinew<-kiwi[,c('id',"scenario",'choice')]
kiwinew1<-merge(kiwinew,product,by='scenario')
kiwinew1<-as.data.frame(kiwinew1[,-13])
kiwinew1<-as.data.frame(sapply(kiwinew1[-1], as.numeric))

calc_shares <- function(coef, x) {
  
  x <- as.matrix(x)
  
  # 1) utility 
  u <-  coef[1] + cbind(x[, 1:3] %*% coef[-1], x[, 4:6] %*% coef[-1], x[, 7:9] %*% coef[-1])
  
  # 2) share = exp / (1 + sum(exp))
  shares <- exp(u) / (1 + rowSums(exp(u))) 
  
  return(shares)
  
}

#the likelihood function
ind.like <- function(coef, x, choice, id) {
  
  # 1) and 2) 
  share <- calc_shares(coef, x)
  
  # 3) likelihood, i.e. the probability of choosing the option in the data
  likelihood = share[, 1] * (choice == 1) + 
    share[, 2] * (choice == 2) + 
    share[, 3] * (choice == 3) + 
    (1 - rowSums(share)) * (choice == 0)
  
  # 4) return product of likelihood within an individual
  indlike <- aggregate(likelihood, by = list(id), FUN = prod)$x
  return(indlike)
  
}

# negative log likelihood with segmentation
latent.like <- function(coef1, coef2, theta, x, choice, id) {
  
  # likelihood for two segments
  l1 <- ind.like(coef1, x, choice, id)
  l2 <- ind.like(coef2, x, choice, id)
  
  # average likelihood
  l <- l1 * theta + l2 * (1 - theta)
  
  # bound l above 0 (practical considerations)
  l <- pmax(l, 1E-200)
  
  # -log likelihood
  return(-sum(log(l)))
  
}

#initial par
coef1 <- c(2,3, 9, -3)
coef2 <- c(2,3, 5, -6)
theta <- 0.2
x <- kiwinew1[,3:11]
choice <- kiwinew1$choice
id <- factor(kiwinew1$id)

# maximum likelihood
latent_seg_est <- optim(par = c(coef1, coef2,theta), fn = function(par) latent.like(par[1:4], par[5:8], par[9], x, choice, id), 
                        method = "L-BFGS-B", lower = c(rep(-Inf, 8), 0), upper = c(rep(Inf, 8), 1), 
                        control = list(pgtol = 1E-12), hessian = TRUE)

# use a function to go over multiple starting point
ml_multi_start <- function(start) {
  est_temp <- optim(par = start, 
                    fn = function(par) latent.like(par[1:4], par[5:8],par[9], x, choice, id),
                    method = "Nelder-Mead", control = list(abstol = 1E-12, reltol = 1E-12))
  return(c(est_temp$par, est_temp$value))
}

# use a set of "reasonable" starting points
set.seed(0)
start_rand <- cbind(runif(15, min = 0, max = 2),runif(15, min = -2, max = 2), runif(15, min = -2, max = 2), runif(15, min = -2, max = 0),
                    runif(15, min = 0, max = 2),runif(15, min =-2, max = 2), runif(15, min = -2, max = 2), runif(15, min = -2, max = 0), 
                    runif(15, min = 0, max = 1))

# use sapply; can also use a for loop
multi_start_check <- sapply(data.frame(t(start_rand)), ml_multi_start)
multi_start_check[,9]



#3[B]segmentation,own and cross-price elasticities
# own price-elasticity and cross-price elasticity
calc_shares_seg <- function(coef_mat, x, segment.share) {
  
  x <- as.matrix(x)
  
  # initialize shares (recommend sticking to row vector x)
  shares <- rep(0, 3)
  
  for (i in 1:length(segment.share)) {
    
    coef <- coef_mat[i, ]
    
    # 1) utility (no intercept for comparison purpose)
    u <- coef[1]+ cbind(x[, 1:3] %*% coef[-1], x[, 4:6] %*% coef[-1], x[, 7:9] %*% coef[-1])
    
    # 2) share = exp / (1 + sum(exp)) for each segment, 
    #   but sum across segments
    shares <- shares + segment.share[i] * exp(u) / (1 + rowSums(exp(u))) 
    
  }
  
  return(shares)
  
}
coef_mat<-t(data.frame(seg1=latent_seg_est$par[1:4],seg2=latent_seg_est$par[5:8]))
segment.share<-c(latent_seg_est$par[9],1-latent_seg_est$par[9])
price_KR<-mean(kiwinew1$price.1)#1.379478
price_KB<-mean(kiwinew1$price.2)#1.375511
price_MB<-mean(kiwinew1[kiwinew1$price.3!=99,"price.3"])#1.345585
char_vec <- t(c(0, 1, price_KR, 1, 1, price_KB, 1, 0,price_MB ))

# perturb prices +0.1(kiwi regular)
s_0 <- calc_shares_seg(coef_mat,char_vec,segment.share)
s_1 <- calc_shares_seg(coef_mat, char_vec + c(0, 0, 0.1, 0, 0, 0,rep(0, 3)),segment.share)

# price elasticities (ask: wait, what are the other things?)
((s_1 - s_0)*price_KB)/(0.1*s_0)

# perturb prices +0.1(kiwi bubble)
s_0 <- calc_shares_seg(coef_mat,char_vec,segment.share)
s_1 <- calc_shares_seg(coef_mat, char_vec + c(0, 0, 0, 0, 0, 0.1,rep(0, 3)),segment.share)

# price elasticities (ask: wait, what are the other things?)
((s_1 - s_0)*price_KB)/(0.1*s_0)

# perturb prices +0.1(mango bubble)
s_0 <- calc_shares_seg(coef_mat,char_vec,segment.share)
s_1 <- calc_shares_seg(coef_mat, char_vec + c(0, 0, 0, 0, 0, 0, 0, 0, 0.1),segment.share)

# price elasticities (ask: wait, what are the other things?)
((s_1 - s_0)*price_KB)/(0.1*s_0)

#Optimal price - before launching kiwi bubbles
mc=0.50
char_vec <- t(c(0, 1, price_KR, 1, 1, 200, 1, 0,1.43 ))
profit <- function(coef, chars, segment.share, col_p, p, mc) {
  
  # note: col_p identifies which column the focal price is at
  chars[col_p] <- p
  
  # calculate shares
  shares_mat <- calc_shares_seg(coef, chars, segment.share)
  
  # calculate profit
  profit <- shares_mat[, 1] * (p - mc)
  
  # return profit
  return(profit)
  
}

profit_max_p <- optim(p = 3, function(p) -profit(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p = 3, p, mc = mc), lower = mc, upper = 10, method = "Brent")

profit_max_p

#Optimal price - after launching kiwi bubbles
mc=0.50
char_vec <- t(c(0, 1, price_KR, 1, 1, price_KB, 1, 0,1.43 ))
profit2 <- function(coef, chars, segment.share, col_p1,col_p2, p1, p2, mc) {
  
  # note: col_p identifies which column the focal price is at
  chars[col_p1] <- p1
  chars[col_p2] <- p2
  
  # calculate shares
  shares_mat <- calc_shares_seg(coef, chars, segment.share)
  
  # calculate profit
  profit2 <- shares_mat[, 1] * (p1 - mc) + shares_mat[, 2] * (p2 - mc)
  
  # return profit
  return(profit2)
  
}
p1<-2
p2<-1
profit_max_p2 <- optim(par = c(p1,p2), 
                       fn=function(par) -profit2(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p1 = 3,col_p2=6, p1=par[1],p2=par[2], mc = mc), 
                       lower = c(mc,mc), upper = c(10,10),control = list(pgtol = 1E-12), hessian = TRUE,method = "L-BFGS-B")
profit_max_p2

#4[C]Testing for consumer behavior and insights for pricing strategies
kiwi$pp<-rep(0,nrow(kiwi))
for (i in 2:nrow(kiwi)) {
  if(kiwi[i,'id']==kiwi[i-1,'id']){
    kiwi[i,'pp']<-kiwi[i-1,'choice']
    if(kiwi[i-1,'choice']==0){
      kiwi[i,'pp']<-kiwi[i-2,'choice']
    }
  }
}
kiwi$pp.1<-rep(0,nrow(kiwi))
kiwi$pp.2<-rep(0,nrow(kiwi))
kiwi$pp.3<-rep(0,nrow(kiwi))
for (i in 1:nrow(kiwi)) {
  if(kiwi[i,'pp']==1){kiwi[i,'pp.1']<-1}
  if(kiwi[i,'pp']==2){kiwi[i,'pp.2']<-1}
  if(kiwi[i,'pp']==3){kiwi[i,'pp.3']<-1}
}

kiwinew2<-data.frame(id=kiwi$id,choice=kiwi$choice,
                     pp.1=kiwi$pp.1, bubble.1=rep(0,nrow(kiwi)),kiwi.1=rep(1,nrow(kiwi)),price.1=kiwi$price_KR,
                     pp.2=kiwi$pp.2, bubble.2=rep(1,nrow(kiwi)),kiwi.2=rep(1,nrow(kiwi)),price.2=kiwi$price_KB,
                     pp.3=kiwi$pp.3, bubble.3=rep(1,nrow(kiwi)),kiwi.3=rep(0,nrow(kiwi)),price.3=kiwi$price_MB)

calc_shares2 <- function(coef, x) {
  
  x <- as.matrix(x)
  
  # 1) utility 
  u <-  coef[1] + cbind(x[, 1:4] %*% coef[2:5], x[, 5:8] %*% coef[2:5], x[, 9:12] %*% coef[2:5])
  
  # 2) share = exp / (1 + sum(exp))
  shares2 <- exp(u) / (1 + rowSums(exp(u))) 
  
  return(shares2)
  
}

#the likelihood function
ind.like2 <- function(coef, x, choice, id) {
  
  # 1) and 2) 
  share <- calc_shares2(coef, x)
  
  # 3) likelihood, i.e. the probability of choosing the option in the data
  likelihood = share[, 1] * (choice == 1) + 
    share[, 2] * (choice == 2) + 
    share[, 3] * (choice == 3) + 
    (1 - rowSums(share)) * (choice == 0)
  
  # 4) return product of likelihood within an individual
  indlike2 <- aggregate(likelihood, by = list(id), FUN = prod)$x
  return(indlike2)
  
}

# negative log likelihood with segmentation
latent.like2 <- function(coef1, coef2, theta, x, choice, id) {
  
  # likelihood for two segments
  l1 <- ind.like2(coef1, x, choice, id)
  l2 <- ind.like2(coef2, x, choice, id)
  
  # average likelihood
  l <- l1 * theta + l2 * (1 - theta)
  
  # bound l above 0 (practical considerations)
  l <- pmax(l, 1E-200)
  
  # -log likelihood
  return(-sum(log(l)))
  
}

#initial par
coef1 <- c(1,1, 1,1,-2)
coef2 <- c(1,1, 1,1, -2)
theta <- 0.2
x <- kiwinew2[,3:14]
choice <- kiwinew2$choice
id <- factor(kiwinew2$id)

# maximum likelihood
latent_seg_est <- optim(par = c(coef1, coef2,theta), fn = function(par) latent.like2(par[1:5], par[6:10], par[11], x, choice, id), 
                        method = "L-BFGS-B", lower = c(rep(-Inf, 10), 0), upper = c(rep(Inf, 10), 1), 
                        control = list(pgtol = 1E-12), hessian = TRUE)


#5.Competition and policy
#strategy actions for Mango
mc=0.50
optkiwi<-profit_max_p2$par
char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0,price_MB ))
profit3 <- function(coef, chars, segment.share, col_p, p, mc) {
  
  # note: col_p identifies which column the focal price is at
  chars[col_p] <- p
  
  # calculate shares
  shares_mat <- calc_shares_seg(coef, chars, segment.share)
  
  # calculate profit
  profit3 <- shares_mat[, 3] * (p - mc)
  
  # return profit
  return(profit3)
  
}
profit_max_p3 <- optim(p = 2, function(p) -profit3(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p = 9, p, mc = mc), lower = mc, upper = 10, method = "Brent")
profit_max_p3
optmango<-profit_max_p3$par
#Strategy action for kiwi if mango set the price at optimal price
char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0,optmango ))
profit_max_p2 <- optim(par = c(p1,p2), 
                       fn=function(par) -profit2(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p1 = 3,col_p2=6, p1=par[1],p2=par[2], mc = mc), 
                       lower = c(mc,mc), upper = c(10,10),control = list(pgtol = 1E-12), hessian = TRUE,method = "L-BFGS-B")
profit_max_p2
optkiwi<-profit_max_p2$par
#Nash Equilibrium
char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0, optmango ))
profit_max_p3 <- optim(p = optmango, function(p) -profit3(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p = 9, p, mc = mc), lower = mc, upper = 10, method = "Brent")
profit_max_p3
optmango<-profit_max_p3$par#0.9030

char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0,optmango ))
profit_max_p2 <- optim(par = c(p1,p2), 
                       fn=function(par) -profit2(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p1 = 3,col_p2=6, p1=par[1],p2=par[2], mc = mc), 
                       lower = c(mc,mc), upper = c(10,10),control = list(pgtol = 1E-12), hessian = TRUE,method = "L-BFGS-B")
profit_max_p2
optkiwi<-profit_max_p2$par#2.6412753 0.9224598

char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0, optmango ))
profit_max_p3 <- optim(p = optmango, function(p) -profit3(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p = 9, p, mc = mc), lower = mc, upper = 10, method = "Brent")
profit_max_p3
optmango<-profit_max_p3$par#0.9026

char_vec <- t(c(0, 1, optkiwi[1], 1, 1, optkiwi[2], 1, 0,optmango ))
profit_max_p2 <- optim(par = c(p1,p2), 
                       fn=function(par) -profit2(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p1 = 3,col_p2=6, p1=par[1],p2=par[2], mc = mc), 
                       lower = c(mc,mc), upper = c(10,10),control = list(pgtol = 1E-12), hessian = TRUE,method = "L-BFGS-B")
profit_max_p2
optkiwi<-profit_max_p2$par#2.6411760 0.9223434

#Collusion and anti_trust policy
profit3 <- function(coef, chars, segment.share, col_p1,col_p2,col_p3,p1, p2, p3, mc) {
  
  # note: col_p identifies which column the focal price is at
  chars[col_p1] <- p1
  chars[col_p2] <- p2
  chars[col_p3] <- p3
  
  # calculate shares
  shares_mat <- calc_shares_seg(coef, chars, segment.share)
  
  # calculate profit
  profit3 <- shares_mat[, 1] * (p1 - mc) + shares_mat[, 2] * (p2 - mc) + shares_mat[, 3] * (p3 - mc)
  
  # return profit
  return(profit3)
  
}
p1<-2
p2<-1
p3<-2
profit_max_p4 <- optim(par = c(p1,p2,p3), 
                       fn=function(par) -profit3(coef = coef_mat, char = char_vec,segment.share=segment.share, col_p1 = 3,col_p2=6,col_p3=9, p1=par[1],p2=par[2],p3=par[3], mc = mc), 
                       lower = c(mc,mc), upper = c(10,10),control = list(pgtol = 1E-12), hessian = TRUE,method = "L-BFGS-B")
profit_max_p4
