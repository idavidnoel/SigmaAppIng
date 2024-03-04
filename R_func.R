center_kernel <- function(x) {
  H <- diag(rep(1, times = ncol(x))) - 1 / ncol(x)
  return(H %*% x %*% H)
}
normalize_kernel <- function(x) {
  D <- diag(1/sqrt(diag(x)))
  return(D %*% x %*% D)
}

scale_kernel_norm <- function(x) {
  return(x / Matrix::norm(x, "F"))
}


K.linear= function(x1){
  #centred by column
  x1<-as.matrix(scale(x1,scale=F))
  KL <- x1%*%(t(x1))
  SV<-mean(diag(KL))
  KL<-KL/SV
  KL
}



K.Polynomial=function(x1, x2=x1, gamma=1/length(x1), b=1, d=3){
  x1<-as.matrix(scale(x1,scale=F))
  PL <- (gamma*(x1%*%t(x1))+b)^d
  SV<-mean(diag(PL))
  PL<-PL/SV
  PL
}

l2norm=function(x){sqrt(sum(x^2))}

K.Gaussian=function(x1,x2=x1, gamma=1/length(x1)){
  #not necessary to center the data
  #x1<-as.matrix(scale(x1,scale=F))
  gauss<-exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol(x2<- t(x2)),
                          Vectorize(function(i, j) l2norm(x1[i,]-x2[,j])^2)))
  colnames(gauss)<-rownames(x1)
  rownames(gauss)<-rownames(x1)
    return(gauss)
}


K.AK1_Final<-function(x1,x2=x1){
  x1<-scale(x1,scale=F)
  x2<-scale(x2,scale=F)
  n1<-nrow(x1)
  n2<-nrow(x2)
  x1tx2<-as.matrix(x1)%*%t(x2)
  norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))
  norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))
  costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)
  costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1
  theta<-acos(costheta)
  normxy<-as.matrix(norm1)%*%t(norm2)
  J = (sin(theta)+(pi-theta)*cos(theta))
  AK1 = 1/pi*normxy*J
  AK1<-AK1/median(AK1)
  colnames(AK1)<-rownames(x2)
  rownames(AK1)<-rownames(x1)
  SV=mean(diag(AK1))
  AK1=AK1/SV #mean diagonal element=1
  return(AK1)
}


MDS_fnc <- function(x1){
  matrix_E <- vegdist(x1, method="bray")
  MDS_mat<-pcoa(matrix_E)
  vec_mds<-as.data.frame(MDS_mat$vectors)
  dat_Mds <- K.linear((vec_mds))
  colnames(dat_Mds)<-rownames(x1)
  rownames(dat_Mds)<-rownames(x1)
  dat_Mds
}


DCA_fnc <- function(x1){
  matrix_E <- vegdist(x1, method="bray")
  dca_mat <- decorana(matrix_E)
  proj_dca <- as.data.frame(dca_mat$rproj)
  dat_DCA <- K.linear((proj_dca))
  colnames(dat_DCA)<-rownames(x1)
  rownames(dat_DCA)<-rownames(x1)
  dat_DCA
}



Euc_fnc <- function(x1){
  matrix_E <- vegdist(x1, method="euclidean")
  matrix_E <- as.matrix(matrix_E)
  #calculate the max diff=similarity of 0.1
  s<-max(matrix_E)
  matrix_E<-matrix(ncol=length(matrix_E[1,]),c(if_else(matrix_E*(1-0.1)/s<1,1-matrix_E*(1-0.1)/s,1)))
  matrix_E
}

PLN_fnc <- function(x1){
  
  x1 <- as.matrix(x1)
  x1 <- t(x1)
  x1 <- as.data.frame(x1)
  x1$cov <- "2"
  cov_tr <- as.data.frame(x1$cov)
  rownames(cov_tr) <- rownames(x1)
  x1<-x1[,!names(x1) %in% c("cov")]
  
  pln_data <- PLNmodels::prepare_data(x1[,1:length(x1)], cov_tr)
  myPLN <- PLN(Abundance~1, data=pln_data)
  PLN_res <- myPLN %>% sigma() 
  PLN_res<-PLN_res/mean(diag(PLN_res))
  PLN_res
}



BC_fnc <- function(x1){
  matrix_E <- vegdist(x1, method="bray")
  matrix_E <- as.matrix(matrix_E)
  matrix_E <- 1-matrix_E
  matrix_E
}


JC_fnc <- function(x1){
  matrix_E <- vegdist(x1, method="jaccard", binary = T)
  matrix_E <- as.matrix(matrix_E)
  matrix_E <- 1-matrix_E
  matrix_E
}

