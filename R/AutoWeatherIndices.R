# This package calculates weather indices from weather variables for using as input for further analysis.
AutoWeatherIndices<-function(x,nw)
{
  imd<-as.matrix(x[,-1])
  yld<-as.numeric(unlist(x[,1]))
  no_of_weeks<-nw
  no_of_w_variables<-round(ncol(imd)/no_of_weeks)
  no_of_years<-as.numeric(nrow(imd))
  wmatrix<- list()
  j<-1
  k<-no_of_weeks
  for (i in 1:no_of_w_variables)
  {
    wmatrix[[i]] <- imd[,j:k]
    j=j+no_of_weeks
    k=k+no_of_weeks
  }

  t<-1:no_of_years;t
  dt<-stats::lm(yld~t)
  adjustedyield<-dt$residuals
  cor_wv<-matrix(nrow=no_of_weeks,ncol=no_of_w_variables,byrow = TRUE)
  varZ1<-matrix(nrow=no_of_years,ncol=no_of_w_variables)
  varZ0<-matrix(nrow=no_of_years,ncol=no_of_w_variables)
  k<-no_of_weeks
  i<-1
  for (j in 1:no_of_w_variables)
  {
    for(i in i:k)
    {
      cor_wv[i]<-stats::cor(adjustedyield,imd[,i])

    }

    i=i+1
    k=k+no_of_weeks
    varZ1[,j]<-(rowSums(wmatrix[[j]]%*%as.vector(cor_wv[,j])))
    varZ0[,j]<-(rowSums(wmatrix[[j]]))
  }
  ##To find the product of the weather variable
  wpmatrix<-list()
  m<-1
  for (p in 1:no_of_w_variables)
  {

    for (q in p+1:no_of_w_variables)
    {
      if(q<=no_of_w_variables)
      {

        wpmatrix[[m]]<-wmatrix[[p]]*wmatrix[[q]]
        m=m+1
      }
    }
  }

  #computing the correlation matrix
  combi<-(gtools::combinations(no_of_w_variables,2))
  tol_combi<-2*nrow(combi)
  corp_wv<-matrix(nrow=no_of_weeks,ncol=nrow(combi),byrow = TRUE)
  varZ11<-matrix(nrow=no_of_years,ncol=nrow(combi))
  varZ10<-matrix(nrow=no_of_years,ncol=nrow(combi))
  k<-no_of_weeks
  i<-1
  l<-1
  for (j in 1:nrow(combi))
  { b=1
  for(i in i:k)
  {
    corp_wv[i]<-stats::cor(adjustedyield,wpmatrix[[j]][,b])
    b=b+1
  }
  i=i+1
  k=k+no_of_weeks
  varZ11[,j]<-(rowSums(wpmatrix[[j]]%*%as.vector(corp_wv[,j])))
  varZ10[,j]<-(rowSums(wpmatrix[[j]]))

  }

  w_indices<-cbind(varZ0,varZ1,varZ10,varZ11)
  correlation_mat_wv<-cor_wv
  correlation_mat_wvp<-corp_wv
  return_list<-list(w_indices,correlation_mat_wv,correlation_mat_wvp)
  return(return_list)
}
