data_prep <- function(data, borrow_hn = NULL){
  # if borrow_hn = NULL, no data borrowing
  # if borrow_hn is a number, then borrow_hn is the bandwidth for data borrowing.
  
    data <- data[unlist(tapply(data$t, data$id, function(x) !duplicated(x))),]  # 1/1000 t is redundant bc rounding
    
    dat <- data[!duplicated(data$id),]
    
    
    data_ls <- c()
    for (i in unique(data$id)){
        data_i <- subset(data,id==i)
        data_ls[[as.character(i)]] <- data_i[,-1]
    }
    
    #---------------------- data with counting process structure -----------------#
    # head(data)
    data_cp <- from_longi_to_cp(longi_data = data, id="id", t="t", X=c("X1", "X2"),
    delta="delta", Y="Y")
    data_cp$start <- round(data_cp$start,3)
    data_cp$stop <- round(data_cp$stop,3)
    
    ### remove those stop<=start, and carefully record delta_star
    idx_rdd <- which(data_cp$stop <= data_cp$start)
    if (length(idx_rdd)>=1) {
        for (k in idx_rdd){
            if (nrow(subset(data_cp, id==data_cp$id[k]))>1){
                ### in case, there is only one obs, then no need to move delta_star
                data_cp$delta_star[k-1] <- data_cp$delta_star[k]
            }
        }
        data_cp <- data_cp[-idx_rdd,]
    }
    
    
    #----------------- data with counting process structure and Xt ---------------#
    temp  <- t(sapply(data_cp$start, fracpoly))
    data_cp <- cbind(data_cp, fracX1=temp[,-1]*data_cp$X1, fracX2=temp[,-1]*data_cp$X2)
    
    #--------------------- proposed with starting value from cox -----------------#
    
    #### step 1 preparation
    # loop over subjects in the list and find out their at-risk subjects
    risk_sets <- list()
    for (i in names(which(lapply(data_ls, function(x) {x[[c("delta")]][1]==1})==TRUE))){
        risk_sets[[i]] =  names(which(lapply(data_ls, function(x) {x[["Y"]][1]>=data_ls[[c(i,"Y")]][1]})==TRUE))
    }
    
    # find out whose risk sets is the subject in? Then find which t's are missing?
    list_for_t_star <- list()
    for (i in names(data_ls)){
        in_whose_risk_set <- names(which(lapply(risk_sets, function(x){sum(x==i)>=1})==TRUE))    # "3" is in risk set of subject "1" and "3"
        
        t_all <- t(unique(subset(data, select=c(t), id%in%in_whose_risk_set)))
        # t_all <- unique(c(sapply(data_ls[c(in_whose_risk_set)], function(x){x[[c("t")]]})))
        
        t_star <- t_all[!t_all%in%data_ls[[c(i, "t")]]]
        
        
        list_for_t_star[[i]] = t_star
        
    }
    
    n_star <- sum(sapply(list_for_t_star, length))
    # do.call(c, lapply(list_for_t_star, length))  // the same
    
    if (is.null(borrow_hn)){
      
      data_cb <- data
      data_cb$group <- 1
      
    } else{
      
      # feed data and t_star to compute X_star
      data_star <- data_ker_C_v3(id=data$id, Y=data$Y, delta=data$delta, X1=data$X1, X2=data$X2,
                                 t=data$t, list_for_t_star=list_for_t_star , h_n=borrow_hn,
                                 id_u= unique(data$id), n_star=n_star)
      
      # same as before, append the X_star, and compute index
      data_star <- as.data.frame(data_star)
      colnames(data_star) <- c("id","delta", "X1", "X2", "t")
      data_cb <- rbind(data[,c("id", "delta", "X1", "X2", "t")],
                       data_star[,c("id", "delta", "X1", "X2", "t")])
      # 1: original measure; 0: missing but later imputed by borrowing info
      data_cb$group <- c(rep(1, nrow(data)), rep(0, nrow(data_star))) 
      
    }
    
    
    data_cb <- data_cb[order(data_cb$id, data_cb$t),]
    index <- gen_index_C_v2(dat_id=dat$id, dat_delta=dat$delta, dat_Y=dat$Y,
    id=data_cb$id, delta=data_cb$delta, t=data_cb$t, X1=data_cb$X1,
    X2=data_cb$X2, group=data_cb$group)
    
    return(list(data_cp=data_cp, data_cb=data_cb, index=index))
}

ker_xtbt <- function(data_cb, index, df, no_marker){
  start.ALL <- rbind(c(1,rep(0,df*no_marker-1)),
  c(-1,rep(0,df*no_marker-1)))
  
  res.ALL <- c()
  C.h <- 1
  for (j in 1:nrow(start.ALL)){
    if( (start.ALL[j,1]) == 1){
      fit <- optim(start.ALL[j,-1], pos_like_two_markers_ker_C, df=df, no_marker=no_marker, X1=data_cb$X1, X2=data_cb$X2, t=data_cb$t, index=index, 
                   h_n=C.h*length(unique(data_cb$id))^(-1/3), control=c(maxit=5000))
      res.ALL <- rbind(res.ALL, c(start=start.ALL[j,], coefs0=1, coefs=fit$par, value=fit$value, counts=fit$counts, convergence0=fit$convergence, message=fit$message))
    }else if ((start.ALL[j,1]) == -1){
      fit <- optim(start.ALL[j,-1], neg_like_two_markers_ker_C, df=df, no_marker=no_marker, X1=data_cb$X1, X2=data_cb$X2, t=data_cb$t, index=index, 
                   h_n=C.h*length(unique(data_cb$id))^(-1/3), control=c(maxit=5000))
      res.ALL <- rbind(res.ALL, c(start=start.ALL[j,], coefs0=-1, coefs=fit$par, value=fit$value, counts=fit$counts, convergence0=fit$convergence, message=fit$message))
    }
    
  }
  
  res.ALL <- data.frame(res.ALL)
  res.ALL <- res.ALL[which.min(res.ALL$value),]
  colnames(res.ALL)[grep("coefs", colnames(res.ALL))] <- c(paste0("X1.", 1:df), paste0("X2.", 1:df))
  
  return(res.ALL)
  
}


from_longi_to_cp <- function(longi_data, id="id", t="v.time", X=c("X1.t", "X2.t"), delta="delta", Y="Y"){
  # the longitudinal data contains id, t, X, delta, Y; sorted by id, t
  # the counting process data contains id, start, time, X, delta_star
  # works for both regularly and irregularly scheduled measurements
  
  longi_data$start <- longi_data[[t]]
  longi_data$stop <- c(longi_data[[t]][-1], 0)
  idx <- cumsum(table(longi_data[,id]))
  longi_data$stop[idx] <- longi_data[[Y]][idx]
  longi_data$delta_star <- 0
  longi_data$delta_star[idx] <- longi_data[[delta]][idx]
  cp_data <- longi_data[,c(id, "start", "stop", "delta_star", X)]
  
  return(cp_data)
  
}

get.ROC.AUC_realdata <- function(data, fit, s, w, 
                                 logscale, offset 
                                 
                                 ){
  ### use IPCW
  
  
  
  ### given those live up to s year,  use the S(X(s);beta(s)) to predict the outcome at the end of s+w yrs.
  ### data must be only have measurements at s, delta is censored at s+w yrs.
  
  ### step 1, get X(s), that is the measurement at landmark time s
  ### 
  sub <- subset(data, Y>s) # exclude subjects who were risk-free at s.
  sub <- subset(sub, t==s)
  sub <- sub[order(sub$id),]  

  sub <- sub[complete.cases(sub),]  # remove those NaN
  
  ### create a new vector to record the survival status at s+w
  sub$delta2 <- sub$delta                                                    
  sub$delta2[sub$Y>(s+w)] <- 0                             # create this new binary outcome; if dead by s+w, then 1.                
  

  
  ### step 2, use SeSp_s_w_c.R to calculate Se at each cutoff, and add up auc empirically.
  gamma <- unlist(fit[grep("X", colnames(fit))])
  ests=list("ker_xtbt"=gamma)
  result <- list()                
  
  
  if(length(table(sub$delta2))==1){
    
    for (i in names(ests)){
      
      result[[i]] <- "no control/case"
      names(result[[i]]) <- "auc"
      
    }
    
  }else{
    
    for (i in names(ests)){
      mtd_name <- i
      
      
      if (i%in%c("ker_xtbt", "cox_xtbt", "pc_xtbt")){
        if (is.na(ests[[i]][1])){
          beta_s <- NA
        } else {
          
          var_nms <- unique(unlist(lapply(strsplit(names(ests[[i]]), "[.]"), function(x) x[1])))
          
          beta_s <- c()
          
          for (var_nm in var_nms){
            var_order <- as.numeric(strsplit(var_nm, "")[[1]][2])
            beta_s <- c(beta_s, 
                        ests[[i]][grep(var_nm, names(ests[[i]]))] %*% fracpoly(s)
            )
          }
          
          
        }
      }
      
      if (is.na(beta_s[1])){
        auc.emp <- NA
      } else if (sum(sub$Y < s+w)==0){  # no event in [s, s+w]
        auc.emp <- NA
      } else {
        
        score <- as.matrix(sub[,grep("X", colnames(sub))])%*%beta_s    # the range of score decide the range of cutoff 
        
        Sp_seq <- Sp_s_w_c_IPCW(Y=sub$Y, delta=sub$delta2, beta=beta_s, X=as.matrix(sub[,grep("X",colnames(sub))]), c=c(min(score)-0.0001,sort(score)), s=s, w=w)
        Se_seq <- Se_s_w_c_IPCW(Y=sub$Y, delta=sub$delta2, beta=beta_s, X=as.matrix(sub[,grep("X",colnames(sub))]), c=c(min(score)-0.0001,sort(score)), s=s, w=w)
        
        if (is.na(Se_seq[1])){
          auc.emp <- NA
        }else{
          idx_axis <- which(!duplicated(Se_seq))
          
          Sp_seq_axix <- Sp_seq[idx_axis]
          Se_seq_axis <- Se_seq[idx_axis]
          
          auc.emp <- as.numeric(diff(Sp_seq_axix)%*%Se_seq_axis[-length(Se_seq_axis)])
        }
        
        
      }
      
      result[[i]] <- c(auc.emp)
      names(result[[i]]) <- c("auc")
    }
    
  }
  
  
  return(result)
  
}

Se_s_w_c_IPCW <- function(Y, delta, beta, X, c, s=0, w=1){
  # Y: observed event time, Y=min(T, C)
  # beta: coefficient;   p
  # X: covariate matrix; n*p
  # c: pre-specified cutoff; vector
  # s: landmark time
  # w: prediction window
  
  score <- as.matrix(X)%*%beta
  
  ### KM estimator for G(t) = Pr(Ci>t)
  cen_surv <- survfit(Surv(Y, 1-delta) ~ 1)
  G_prob <- sapply(Y , function(x) summary(cen_surv, time=x)$surv)
  
  if (length(which(G_prob==0))>=1){
    idx <- which(G_prob==0)
    
    if (sum(delta[idx] == 0)>=1) {          # if for the censor variable, the last subject had event, then must 1-delta=1
      G_prob[idx] = 0.001           # so this change has no effect.
      
    }else{
      cat("ERROR in Gprob \n")
      
    }
  }
  
  IPCW <- delta/G_prob
  
  
  ### re-weight contribution of uncensored subjects using IPCW
  # Se = sum(1*(score > c)*(Y >= s)*(Y <= s+w)*IPCW) / sum(1*((Y >= s) & (Y <= s+w))*IPCW)
  Se = sapply(c, function(x) sum(1*(score > x)*(Y >= s)*(Y <= s+w)*IPCW) / sum(1*((Y >= s) & (Y <= s+w))*IPCW))
  
  
  # Sp = sum(1*(score <= c)*(Y > s+w)*IPCW) / sum(1*(Y > s+w)*IPCW)   # This is wrong
  # Sp = sum(1*(score <= c)*(Y > s+w)) / sum(1*(Y > s+w))               # This is correct, Se Uno, Cai 2007 JASA
  
  
  # return(list(Se=Se, Sp=Sp))
  return(Se)
}

Sp_s_w_c_IPCW <- function(Y, delta, beta, X, c, s=0, w=1){
  # Y: observed event time, Y=min(T, C)
  # beta: coefficient;   p
  # X: covariate matrix; n*p
  # c: pre-specified cutoff; vector
  # s: landmark time
  # w: prediction window
  
  score <- as.matrix(X)%*%beta
  
  # Sp = sum(1*(score <= c)*(Y > s+w)) / sum(1*(Y > s+w))               # This is correct, See Uno, Cai 2007 JASA
  Sp <- sapply(c, function(x) sum(1*(score <= x)*(Y > s+w)) / sum(1*(Y > s+w))      )
  
  return(Sp)
}


get.D.s_realdata <- function(data, fit, s, logscale, offset,
                             tau.rmst
                             ){
  ### step 1, get X(s), that is the measurement at landmark time s
  sub <- subset(data, Y>s) # exclude subjects who were risk-free at s.
  sub <- subset(sub, t==s)
  sub <- sub[order(sub$id),]  
  sub <- sub[complete.cases(sub),]  # remove those NaN
  
  

  
  ### Step 2, use package to calculate RMST.
  gamma <- unlist(fit[grep("X", colnames(fit))])
  ests=list("ker_xtbt"=gamma)
  result <- list()
  
  for (i in names(ests)){
    
    mtd_name <- i
    
    
    if (i%in%c("ker_xtbt", "cox_xtbt", "pc_xtbt")){
      if (is.na(ests[[i]][1])){
        beta_s <- NA
      } else {
        
        var_nms <- unique(unlist(lapply(strsplit(names(ests[[i]]), "[.]"), function(x) x[1])))
        
        beta_s <- c()
        
        for (var_nm in var_nms){
          var_order <- as.numeric(strsplit(var_nm, "")[[1]][2])
          beta_s <- c(beta_s, 
                      ests[[i]][grep(var_nm, names(ests[[i]]))] %*% fracpoly(s)
          )
        }
        
        
      }
    }
    
    
    
    if(is.na(beta_s[1])){
      D_rmst <- NA
    }else{
      score_s <- as.matrix(sub[,grep("X", colnames(sub))])%*%beta_s
      
      group <- ifelse(score_s>median(score_s), 1, 0)
      rmst <- rmst2(time=sub$Y, status=sub$delta, arm=group, tau=tau.rmst)
      D_rmst <- -rmst$unadjusted.result[1,1]
    }
    
    result[[i]] <- D_rmst 
    names(result[[i]]) <- c("d")
  }
  
  return(result)
}

get.score <- function(data, fit, s, df){
# df: df of fractional polynomials
  
  
  sub <- subset(data, t==s)
  s_frac <- fracpoly(s)
  
  gamma <- unlist(fit[grep("X", colnames(fit))])
  
  as.matrix(sub[,grep("X", colnames(sub))])%*%(matrix(gamma[c(grep("X1",names(gamma)),
                                                                     grep("X2",names(gamma)))], 
                                                      ncol=df.fracpoly, 
                                                      byrow=TRUE)%*%s_frac)
  
}
  