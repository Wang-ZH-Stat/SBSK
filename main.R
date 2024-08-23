
library(MASS)
library(GeneralizedHyperbolic)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
source("functions.R")

# Explore feasible scope of T2
avg_t_T2 <- c()
for(b_r in c(1:20)){
  avg_t_T2[b_r] <- get_b_T2(ARL = 5000, B0 = 20, range = b_r + 30)[2]
}

# Explore feasible scope of GLR
avg_t_GLR <- c()
for(b_r in c(1:20)){
  avg_t_GLR[b_r] <- get_b_GLR(ARL = 5000, B0 = 20, range = b_r + 115)[2]
}

# Thresholds for three methods
b_scan <- get_b_scan(ARL = 5000, B0 = 20)
b_T2 <- get_b_T2(ARL = 5000, B0 = 20, range = c(31, 50))
b_GLR <- get_b_GLR(ARL = 5000, B0 = 20, range = c(116, 135))

# Experiments
# Experiment I
exper_res <- function(case){
  set.seed(9)
  count_list <- c(0, 0, 0)
  DD_scan <- c()
  DD_T2 <- c()
  DD_GLR <- c()
  for(i in c(1:5000)){
    print(paste('i:', i))
    if(case == 1){
      # mean shift
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = diag(10))
    }
    else if(case == 2){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      Sigma_2 <- diag(10)
      for(j in c(1:5)){
        Sigma_2[j,j] <- 2
      }
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = Sigma_2)
    }
    else if(case == 3){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = 2*diag(10))
    }
    else if(case == 4){
      # Gaussian mixture
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- 0.3*mvrnorm(100, mu = rep(0,10), Sigma = 0.1*diag(10)) + 
        0.7*mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
    }
    else if(case == 5){
      # Gaussian to Laplace
      X1 <- as.matrix(rnorm(100))
      X2 <- as.matrix(rskewlap(100, mu = 0, alpha = 1, beta = 1))
    }
    X <- rbind(X1,X2)
    
    if(count_list[1] < 500){
      t_scan <- dete_scan(X, N = 5, B0 = 20, b = b_scan, sigma = median(dist(X)))[1]
      if(t_scan <= 150){
        count_list[1] <- count_list[1] + 1
        DD_scan[count_list[1]] <- t_scan
        print(paste('count_list_scan:',count_list[1]))
      }
    }
    
    if(count_list[2] < 500){
      t_T2 <- dete_T2(X, B0 = 20, b = b_T2)[1]
      if(t_T2 <= 150){
        count_list[2] <- count_list[2] + 1
        DD_T2[count_list[2]] <- t_T2
        print(paste('count_list_T2:',count_list[2]))
      }
    }
    
    if(count_list[3] < 500){
      t_GLR <- dete_GLR(X, B0 = 20, b = b_GLR)[1]
      if(t_GLR <= 150){
        count_list[3] <- count_list[3] + 1
        DD_GLR[count_list[3]] <- t_GLR
        print(paste('count_list_GLR:',count_list[3]))
      }
    }
    
    if(sum(count_list) == 1500){
      break
    }
  }
  
  return(list(count_list, DD_scan, DD_T2, DD_GLR))
}

exp_res_1 <- exper_res(case = 1)
exp_res_2 <- exper_res(case = 2)
exp_res_3 <- exper_res(case = 3)
exp_res_4 <- exper_res(case = 4)
exp_res_5 <- exper_res(case = 5)

# Show results
round(c(exp_res_1[[1]], mean(exp_res_1[[2]]), mean(exp_res_1[[3]]), mean(exp_res_1[[4]])),2)
round(c(exp_res_2[[1]], mean(exp_res_2[[2]]), mean(exp_res_2[[3]]), mean(exp_res_2[[4]])),2)
round(c(exp_res_3[[1]], mean(exp_res_3[[2]]), mean(exp_res_3[[3]]), mean(exp_res_3[[4]])),2)
round(c(exp_res_4[[1]], mean(exp_res_4[[2]]), mean(exp_res_4[[3]]), mean(exp_res_4[[4]])),2)
round(c(exp_res_5[[1]], mean(exp_res_5[[2]]), mean(exp_res_5[[3]]), mean(exp_res_5[[4]])),2)
round(c(exp_res_1[[1]], sd(exp_res_1[[2]]), sd(exp_res_1[[3]]), sd(exp_res_1[[4]])),2)
round(c(exp_res_2[[1]], sd(exp_res_2[[2]]), sd(exp_res_2[[3]]), sd(exp_res_2[[4]])),2)
round(c(exp_res_3[[1]], sd(exp_res_3[[2]]), sd(exp_res_3[[3]]), sd(exp_res_3[[4]])),2)
round(c(exp_res_4[[1]], sd(exp_res_4[[2]]), sd(exp_res_4[[3]]), sd(exp_res_4[[4]])),2)
round(c(exp_res_5[[1]], sd(exp_res_5[[2]]), sd(exp_res_5[[3]]), sd(exp_res_5[[4]])),2)

# Plot
plot_com_method <- data.frame(Case = c(rep("Case 1", 1500), rep("Case 2", 1500), 
                                       rep("Case 3", 1500), rep("Case 4", 1000), 
                                       rep("Case 5", 500)), 
                              Statistic = c(rep(c(rep("Scan B", 500), rep("T2", 500), rep("GLR", 500)), 3), 
                                            rep("Scan B", 500), rep("GLR", 500), rep("Scan B", 500)), 
                              Delay = c(exp_res_1[[2]]-100, exp_res_1[[3]]-100, exp_res_1[[4]]-100, 
                                        exp_res_2[[2]]-100, exp_res_2[[3]]-100, exp_res_2[[4]]-100, 
                                        exp_res_3[[2]]-100, exp_res_3[[3]]-100, exp_res_3[[4]]-100, 
                                        exp_res_4[[2]]-100, exp_res_4[[4]]-100, exp_res_5[[2]]-100))

plot_com_method <- data.frame(Case = c(rep("Case 1", 1500), rep("Case 2", 1500), 
                                       rep("Case 3", 1500), rep("Case 4", 1500), 
                                       rep("Case 5", 1500)), 
                              Statistic = rep(c(rep("Scan B", 500), rep("T2", 500), rep("GLR", 500)), 5), 
                              Delay = c(exp_res_1[[2]]-100, exp_res_1[[3]]-100, exp_res_1[[4]]-100, 
                                        exp_res_2[[2]]-100, exp_res_2[[3]]-100, exp_res_2[[4]]-100, 
                                        exp_res_3[[2]]-100, exp_res_3[[3]]-100, exp_res_3[[4]]-100, 
                                        exp_res_4[[2]]-100, rep(50, 500), exp_res_4[[4]]-100,
                                        exp_res_5[[2]]-100, rep(50, 500), rep(50, 500)))

ggplot(data=plot_com_method, aes(x=Case, y=Delay, fill=Statistic)) + geom_boxplot() + 
  theme_bw() + theme(legend.position = c(0.1,0.85), axis.title.x=element_text(size=13), 
                     axis.title.y=element_text(size=13), title=element_text(size=12, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_fill_brewer(palette = "Set3")


# Experiment II
exper_com_N <- function(case, N_range, iter = 100){
  set.seed(9)
  num_N <- length(N_range)
  count_list <- rep(0, num_N)
  EDD_list <- matrix(rep(0.0, num_N*iter), nrow = num_N)
  
  for(i in c(1:5000)){
    print(paste('i:', i))
    if(case == 1){
      # mean shift
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = diag(10))
    }
    else if(case == 2){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      Sigma_2 <- diag(10)
      for(j in c(1:5)){
        Sigma_2[j,j] <- 2
      }
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = Sigma_2)
    }
    else if(case == 3){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = 2*diag(10))
    }
    else if(case == 4){
      # Gaussian mixture
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- 0.3*mvrnorm(100, mu = rep(0,10), Sigma = 0.1*diag(10)) + 
        0.7*mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
    }
    else if(case == 5){
      # Gaussian to Laplace
      X1 <- as.matrix(rnorm(100))
      X2 <- as.matrix(rskewlap(100, mu = 0, alpha = 1, beta = 1))
    }
    X <- rbind(X1,X2)
    
    for(j in c(1:num_N)){
      print(paste('N:',N_range[j]))
      if(count_list[j] < iter){
        t_scan <- dete_scan(X, N = N_range[j], B0 = 20, b = b_scan, sigma = median(dist(X)))[1]
        if(t_scan <= 150){
          count_list[j] <- count_list[j] + 1
          EDD_list[j, count_list[j]] <- t_scan
          print(paste('count_list_', j, ':',count_list[j]))
        }
      }
    }
    
    if(sum(count_list) == 100*num_N){
      break
    }
  }
  
  return(EDD_list)
}

N_range <- c(1,5,9,13,17,21)
com_N_res_1 <- exper_com_N(case = 1, N_range = N_range)
com_N_res_2 <- exper_com_N(case = 2, N_range = N_range)
com_N_res_3 <- exper_com_N(case = 3, N_range = N_range)
com_N_res_4 <- exper_com_N(case = 4, N_range = N_range)
com_N_res_5 <- exper_com_N(case = 5, N_range = N_range)

round(apply(com_N_res_1, 1, mean), 2)
round(apply(com_N_res_2, 1, mean), 2)
round(apply(com_N_res_3, 1, mean), 2)
round(apply(com_N_res_4, 1, mean), 2)
round(apply(com_N_res_5, 1, mean), 2)
round(apply(com_N_res_1, 1, sd), 2)
round(apply(com_N_res_2, 1, sd), 2)
round(apply(com_N_res_3, 1, sd), 2)
round(apply(com_N_res_4, 1, sd), 2)
round(apply(com_N_res_5, 1, sd), 2)

plot_com_N <- data.frame(rep(c(N_range),3), c(rep("Case 1", 6), rep("Case 2", 6), rep("Case 3", 6)),
                         c(apply(com_N_res_1, 1, mean)-100, 
                           apply(com_N_res_2, 1, mean)-100, 
                           apply(com_N_res_3, 1, mean)-100), 
                         c(apply(com_N_res_1, 1, sd), 
                           apply(com_N_res_2, 1, sd), 
                           apply(com_N_res_3, 1, sd)))
names(plot_com_N) <- c("N", "Case", "EDD", "sd")
g1 <- ggplot(data = plot_com_N, aes(x=N, y=EDD, group = Case, color= Case, shape= Case)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + geom_errorbar(aes(ymin=EDD-sd,ymax=EDD+sd), width = 0.5) +
  xlab("N")+ ylab("Expected Detection Delay") + 
  labs(title = "Panel A: Different number of reference blocks") + 
  theme_bw() + theme(legend.position = c(0.8,0.85), axis.title.x=element_text(size=13), 
                     axis.title.y=element_text(size=13), title=element_text(size=12, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(1,21),breaks = seq(1,21,4))

# Experiment III
exper_com_B0 <- function(case, B0_range, iter = 100){
  set.seed(9)
  num_B0 <- length(B0_range)
  count_list <- rep(0, num_B0)
  EDD_list <- matrix(rep(0, num_B0*iter), nrow = num_B0)
  
  for(i in c(1:5000)){
    print(paste('i:', i))
    if(case == 1){
      # mean shift
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = diag(10))
    }
    else if(case == 2){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      Sigma_2 <- diag(10)
      for(j in c(1:5)){
        Sigma_2[j,j] <- 2
      }
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = Sigma_2)
    }
    else if(case == 3){
      # covariance change
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(1,10), Sigma = 2*diag(10))
    }
    else if(case == 4){
      # Gaussian mixture
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- 0.3*mvrnorm(100, mu = rep(0,10), Sigma = 0.1*diag(10)) + 
        0.7*mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
    }
    else if(case == 5){
      # Gaussian to Laplace
      X1 <- as.matrix(rnorm(100))
      X2 <- as.matrix(rskewlap(100, mu = 0, alpha = 1, beta = 1))
    }
    X <- rbind(X1,X2)
    
    for(j in c(1:num_B0)){
      print(paste('B0:',B0_range[j]))
      b_scan_t <- get_b_scan(ARL = 5000, B0 = B0_range[j])
      if(count_list[j] < iter){
        t_scan <- dete_scan(X, N = 5, B0 = B0_range[j], b = b_scan_t, sigma = median(dist(X)))[1]
        if(t_scan <= 150){
          count_list[j] <- count_list[j] + 1
          EDD_list[j, count_list[j]] <- t_scan
          print(paste('count_list_', j, ':',count_list[j]))
        }
      }
    }
    
    if(sum(count_list) == iter*num_B0){
      break
    }
  }
  
  return(EDD_list)
}

B0_range <- c(5,10,15,20,25,30)
com_B0_res_1 <- exper_com_B0(case = 1, B0_range = B0_range)
com_B0_res_2 <- exper_com_B0(case = 2, B0_range = B0_range)
com_B0_res_3 <- exper_com_B0(case = 3, B0_range = B0_range)
com_B0_res_4 <- exper_com_B0(case = 4, B0_range = B0_range)
com_B0_res_5 <- exper_com_B0(case = 5, B0_range = B0_range)

round(apply(com_B0_res_1, 1, mean), 2)
round(apply(com_B0_res_2, 1, mean), 2)
round(apply(com_B0_res_3, 1, mean), 2)
round(apply(com_B0_res_4, 1, mean), 2)
round(apply(com_B0_res_5, 1, mean), 2)
round(apply(com_B0_res_1, 1, sd), 2)
round(apply(com_B0_res_2, 1, sd), 2)
round(apply(com_B0_res_3, 1, sd), 2)
round(apply(com_B0_res_4, 1, sd), 2)
round(apply(com_B0_res_5, 1, sd), 2)

plot_com_B0 <- data.frame(rep(c(B0_range),3), c(rep("Case 1", 6), rep("Case 2", 6), rep("Case 3", 6)),
                          c(apply(com_B0_res_1, 1, mean)-100, 
                            apply(com_B0_res_2, 1, mean)-100, 
                            apply(com_B0_res_3, 1, mean)-100), 
                          c(apply(com_B0_res_1, 1, sd), 
                            apply(com_B0_res_2, 1, sd), 
                            apply(com_B0_res_3, 1, sd)))
names(plot_com_B0) <- c("B0", "Case", "EDD", "sd")
g2 <- ggplot(data = plot_com_B0, aes(x=B0, y=EDD, group = Case, color= Case, shape= Case)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + geom_errorbar(aes(ymin=EDD-sd,ymax=EDD+sd), width = 0.5) + 
  xlab(expression(B[0]))+ ylab("Expected Detection Delay") + 
  labs(title = "Panel B: Different block sizes") + 
  theme_bw() + theme(legend.position = c(0.2,0.85), axis.title.x=element_text(size=12), 
                     axis.title.y=element_text(size=13), title=element_text(size=13, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(5,30),breaks = seq(5,30,5))

# Plot
grid.arrange(g1, g2, ncol=2)


# Experiment IV
exper_com_sigma_mu <- function(sigma_range, mu_range, iter = 100, kernel = "Gaussian"){
  set.seed(9)
  num_sigma <- length(sigma_range)
  num_mu <- length(mu_range)
  count_list <- matrix(rep(0, num_sigma*num_mu), nrow = num_sigma)
  EDD_list <- matrix(rep(0, num_sigma*num_mu), nrow = num_sigma)
  for(l in c(1:(10*iter))){
    print(paste('iter:', l))
    for(i in c(1:num_mu)){
      print(paste('i:', i))
      X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
      X2 <- mvrnorm(100, mu = rep(mu_range[i],10), Sigma = diag(10))
      X <- rbind(X1,X2)
      for(j in c(1:num_sigma)){
        print(paste('j:', j))
        if(count_list[j, i] < iter){
          t_scan <- dete_scan(X, N = 5, B0 = 20, b = b_scan, sigma = sigma_range[j]*median(dist(X)), kernel = kernel)[1]
          if(t_scan <= 150){
            count_list[j, i] <- count_list[j, i] + 1
            EDD_list[j, i] <- EDD_list[j, i] + t_scan / iter
            print(paste('count_list_', j,i, ':',count_list[j,i]))
          }
        }
      }
    }
    
    if(sum(count_list) == iter*num_mu*num_sigma){
      break
    }
  }
  
  return(EDD_list)
}

sigma_range <- c(0.5, 1, 10, 100)
mu_range <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

# Gaussian kernel
com_sigma_mu_gau <- exper_com_sigma_mu(sigma_range = sigma_range, mu_range = mu_range)
plot_sigma_mu_gau <- data.frame(sigma = c(rep("0.5 median", 10), rep("1.0 median", 10), 
                                          rep("10 median", 10), rep("100 median", 10)), 
                                mu = rep(mu_range, 4), 
                                EDD = c(com_sigma_mu_gau[1,]-100, com_sigma_mu_gau[2,]-100, 
                                        com_sigma_mu_gau[3,]-100, com_sigma_mu_gau[4,]-100))

g3 <- ggplot(data = plot_sigma_mu_gau, aes(x=mu, y=EDD, group = sigma, color= sigma, shape= sigma)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + 
  xlab(expression(mu))+ ylab("Expected Detection Delay") + 
  labs(title = "Panel A: Gaussian kernel", fill=expression(sigma)) + 
  theme_bw() + theme(legend.position = c(0.8,0.75), axis.title.x=element_text(size=12), 
                     axis.title.y=element_text(size=13), title=element_text(size=13, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(0.2,2.0),breaks = seq(0.2,2,0.2)) + 
  scale_y_continuous(limits = c(3,21),breaks = seq(3,21,3))

# Laplacian kernel
com_sigma_mu_lap <- exper_com_sigma_mu(sigma_range = sigma_range, mu_range = mu_range, kernel = "Laplacian")
plot_sigma_mu_lap <- data.frame(sigma = c(rep("0.5 median", 10), rep("1.0 median", 10), 
                                          rep("10 median", 10), rep("100 median", 10)), 
                                mu = rep(mu_range, 4), 
                                EDD = c(com_sigma_mu_lap[1,]-100, com_sigma_mu_lap[2,]-100, 
                                        com_sigma_mu_lap[3,]-100, com_sigma_mu_lap[4,]-100))
g4 <- ggplot(data = plot_sigma_mu_lap, aes(x=mu, y=EDD, group = sigma, color= sigma, shape= sigma)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + 
  xlab(expression(mu))+ ylab("Expected Detection Delay") + 
  labs(title = "Panel B: Laplacian kernel", fill=expression(sigma)) + 
  theme_bw() + theme(legend.position = c(0.8,0.75), axis.title.x=element_text(size=12), 
                     axis.title.y=element_text(size=13), title=element_text(size=13, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(0.2,2.0),breaks = seq(0.2,2,0.2))+ 
  scale_y_continuous(limits = c(3,21),breaks = seq(3,21,3))

# Plot
grid.arrange(g3, g4, ncol=2)

# Experiment V
# Improvement
exper_com_im_ARL <- function(ARL_range, iter = 100){
  num_ARL <- length(ARL_range)
  count_list <- matrix(rep(0, num_ARL*2), nrow = 2)
  EDD_list_1 <- matrix(rep(0, num_ARL*iter), nrow = num_ARL)
  EDD_list_2 <- matrix(rep(0, num_ARL*iter), nrow = num_ARL)
  set.seed(1)
  for(i in c(1:5000)){
    print(paste("i:",i))
    X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
    X1 <- mvrnorm(100, mu = rep(1,10), Sigma = diag(10))
    X <- rbind(X1,X2)
    for(j in c(1:num_ARL)){
      print(paste("j:",j))
      b_scan_t <- get_b_scan(ARL = ARL_range[j], B0=20)
      if(count_list[1,j] < iter){
        t_scan <- dete_scan(X, N = 5, B0 = 20, b = b_scan_t, sigma = median(dist(X)))[1]
        if(t_scan <= 150){
          count_list[1,j] <- count_list[1,j] + 1
          EDD_list_1[j, count_list[1,j]] <- t_scan
          print(paste('count_list_', 1,j, ':',count_list[1]))
        }
      }
      
      if(count_list[2,j] < iter){
        t_scan <- dete_scan(X, N = 5, B0 = 20, b = b_scan_t, sigma = median(dist(X)), improve = TRUE)[1]
        if(t_scan <= 150){
          count_list[2,j] <- count_list[2,j] + 1
          EDD_list_2[j, count_list[2,j]] <- t_scan
          print(paste('count_list_', 2,j, ':',count_list[2]))
        }
      }
    }
    
    if(sum(count_list) == 2*iter*num_ARL){
      break
    }
  }
  
  return(list(EDD_list_1, EDD_list_2))
}

B0_range <- c(5,10,15,20,25,30)
ARL_range <- 1000*c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Explore the value of thresholds
b_scan_list <- matrix(rep(0, 6*10), nrow = 6)
for(i in c(1:6)){
  for(j in c(1:10)){
    b_scan_list[i,j] <- get_b_scan(ARL_range[j], B0_range[i])
  }
}
plot_b_scan <- data.frame(ARL = rep(ARL_range,6), 
                          B0 = as.factor(c(rep(5, 10), rep(10, 10), rep(15, 10), 
                                           rep(20, 10), rep(25, 10), rep(30, 10))), 
                          b_scan = c(b_scan_list[1,], b_scan_list[2,], 
                                     b_scan_list[3,], b_scan_list[4,], 
                                     b_scan_list[5,], b_scan_list[6,]))
g5 <- ggplot(data = plot_b_scan, aes(x=ARL, y=b_scan, group = B0, color= B0, shape= B0)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + 
  xlab("Average Run Length")+ ylab("Threshold b") + 
  labs(title = "Panel A: Theoretical thresholds", group="B0") + 
  theme_bw() + theme(legend.position = c(0.8,0.25), axis.title.x=element_text(size=12), 
                     axis.title.y=element_text(size=13), title=element_text(size=13, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(1000,10000),breaks = seq(1000,10000,1000), 
                     labels = c("1k","2k", "3k", "4k","5k", "6k", "7k", "8k", "9k", "10k"))

com_im_ARL_res <- exper_com_im_ARL(ARL_range, iter = 100)
non_impr_mean <- apply(com_im_ARL_res[[1]], 1, mean)
impr_mean <- apply(com_im_ARL_res[[2]], 1, mean)

plot_com_ARL <- data.frame(ARL = rep(ARL_range,2), 
                           Subsampling = c(rep("Completely random", 10), rep("Incompletely random", 10)), 
                                        EDD = c(non_impr_mean-100, impr_mean-100))
g6 <- ggplot(data = plot_com_ARL, aes(x=ARL, y=EDD, group = Subsampling, color= Subsampling, shape= Subsampling)) +
  geom_point(size=2) + geom_line(lwd = 1.3) + 
  xlab("Average Run Length")+ ylab("Expected Detection Delay") + 
  labs(title = "Panel B: Comparison of different subsampling methods") + 
  theme_bw() + theme(legend.position = c(0.8,0.25), axis.title.x=element_text(size=12), 
                     axis.title.y=element_text(size=13), title=element_text(size=13, face="italic",hjust=0.5), 
                     plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(1000,10000),breaks = seq(1000,10000,1000), 
                     labels = c("1k","2k", "3k", "4k","5k", "6k", "7k", "8k", "9k", "10k"))

# Plot
grid.arrange(g5, g6, ncol=2)

# Others
X1 <- mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
X2 <- 0.3*mvrnorm(100, mu = rep(0,10), Sigma = 0.1*diag(10)) + 
  0.7*mvrnorm(100, mu = rep(0,10), Sigma = diag(10))
X <- rbind(X1,X2)
dete_scan(X, N = 5, B0 = 20, b = b_scan, sigma = median(dist(X)))
dete_scan(X, N = 5, B0 = 20, b = b_scan, sigma = median(dist(X)), improve = TRUE)
dete_T2(X, B0 = 20, b = b_T2)
dete_GLR(X, B0 = 20, b = b_GLR)



