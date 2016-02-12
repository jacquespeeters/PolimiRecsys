library(readr)
library(Matrix)
library(dplyr)
library(rARPACK)

# Mes fonctions
##################################
# Fonctions qui permettent le calcul MAP
apk <- function(k, actual, predicted){
  score <- 0.0
  cnt <- 0.0
  for (i in 1:min(k,length(predicted)))
  {
    if (predicted[i] %in% actual && !(predicted[i] %in% predicted[0:(i-1)]))
    {
      cnt <- cnt + 1
      score <- score + cnt/i 
    }
  }
  score <- score / min(length(actual), k)
  score
}

mapk <- function (k, actual, predicted){
  if( length(actual)==0 || length(predicted)==0 ) 
  {
    return(0.0)
  }
  
  scores <- rep(0, length(actual))
  for (i in 1:length(scores))
  {
    scores[i] <- apk(k, actual[[i]], predicted[[i]])
  }
  score <- mean(scores)
  score
}

# Partie URMsvd
URM = function(train){
  cat("Calcul de l'URM")
  sparseTrain = sparseMatrix(unlist(train[,1]), unlist(train[,2]), x = unlist(train[,3]))
  return(sparseTrain)
}

maSvd=function(sparseTrain, eigen){
  cat("\n\nCalcul de la svd\n")
  svd = svds(sparseTrain, eigen, nu = eigen, nv = eigen)
  return (svd)
}

URMsvd = function(train, iterUser, Ndim){
  sparseTrain= URM(train)
  
  sparseTrainR1=sparseTrain
  sparseTrainR1@x[]=1
  
  svd=maSvd(sparseTrain,Ndim)
  
  score = as.matrix(svd$u[iterUser,1:Ndim] %*% Diagonal(x=svd$d[1:Ndim]) %*% t(svd$v[,1:Ndim]))
  score = score -1000*as.matrix(sparseTrainR1[iterUser,])
  
  submission=matrix(data = NA, nrow = length(iterUser), ncol = 2)
  colnames(submission) <- c("userId","testItems")
  submission = as.data.frame(submission)
  submission$userId= iterUser
  submission$testItems = do.call(rbind, lapply(1:length(iterUser),FUN=function(i){
    # cat(i," ")
    res = list(order(score[i,], decreasing = TRUE)[1:5]) #Vérifier si ressort bien les noms des colonnes et non les indices
    return(res)
  }))
  return (submission)
}

URMsvdBest = function(train, iterUser, iterDim){
  sparseTrain= URM(train)
  
  sparseTrainR1=sparseTrain
  sparseTrainR1@x[]=1
  
  submission=matrix(data = NA, nrow = length(iterUser), ncol = 2)
  colnames(submission) <- c("userId","testItems")
  submission = as.data.frame(submission)
  submission$userId= iterUser
  
  svd=maSvd(sparseTrain,max(iterDim))
  
  for (i in iterDim){
    score = as.matrix(svd$u[iterUser,1:i] %*% Diagonal(x=svd$d[1:i]) %*% t(svd$v[,1:i]))
    score = score -1000*as.matrix(sparseTrainR1[iterUser,])
    k = do.call(rbind, lapply(1:length(iterUser),FUN=function(i){
      # cat(i," ")
      res = list(order(score[i,], decreasing = TRUE)[1:5]) #Vérifier si ressort bien les noms des colonnes et non les indices
      return(res)
    }))
    submission = cbind(submission, k)
  }
  return (submission)
}

URMsvdBest2 = function(train, iterDim, nUsers){
  train$index = 1:nrow(train)
  usersLocal = train %>% filter(rating>7) %>% group_by(userId) %>% summarise(nb = n()) %>% filter(nb > 1)
  usersLocal = sample(usersLocal$userId, nUsers)
  
  testLocal = train %>% group_by(itemId) %>% mutate( pop = n()) %>%filter(rating>7, userId %in% usersLocal) %>% group_by(userId) %>% arrange(pop)%>% slice(1:(n()/2))    
  trainLocal = train %>% filter(!index %in% testLocal$index)
  
  usersLocal =  sort(usersLocal)
  actual = testLocal %>% group_by(userId) %>% summarise(actual= list(itemId) )
  
  submission = URMsvdBest(trainLocal, usersLocal, iterDim)
  submission$actual = actual$actual
  return(submission)
}

URMsvdBest3 = function( train, iterDim, n, nUsers){
  sub = c()
  for ( i in 1:n){
    t = proc.time()
    cat("Etape", i," du bootstrap BestSVd \n")
    sub = rbind(sub,URMsvdBest2(train, iterDim, nUsers))
    cat( round((proc.time() - t)[3]/60, 1), "mins\n")
  }
  res = c()
  for ( i in 1:length(iterDim)){
    res=c(res, mapk(10, sub[,(i+2)], sub$actual))
  }
  names(res) = iterDim
  return(res)
}

# Partie ContentBased
contentBased = function(train, iterUser){
  sparseTrain = sparseMatrix(unlist(train[,1]), unlist(train[,2]), x = unlist(train[,3]))

  # Création de sparseIcm
  icm=read_csv("icm.csv")
  icm[,3]=1
  sparseIcm = sparseMatrix(unlist(icm[,2]), unlist(icm[,1]), x = unlist(icm[,3]))
  
  # Il faut que les deux sparse soient de même dim...
  # sparseIcm=t(sparseIcm) #La mettre à l'endroit
  vectCol = 1:min(sparseTrain@Dim[2],sparseIcm@Dim[2])
  sparseTrain=sparseTrain[,vectCol]
  sparseTrainR1 = sparseTrain
  sparseTrainR1@x[] = 1
  sparseIcm=sparseIcm[,vectCol]
  rm(vectCol, icm)
  
#   svd=maSvd(sparseTrain,Ndim)
#   sparseTrainSvd = as.matrix(svd$u[iterUser,1:Ndim] %*% Diagonal(x=svd$d[1:Ndim]) %*% t(svd$v[,1:Ndim]))
  
  idf = as.vector( log(ncol(sparseIcm)  / rowSums(sparseIcm)) )
  sparseIcm = (sparseIcm) * idf
  rm(idf)
  
  sparseTrain = sparseTrain[iterUser,]
  
  # Normaliser sparseICM
  norm = sqrt(colSums(sparseIcm^2))
  norm[which(norm == 0)]=1
  sparseIcm = t(t(sparseIcm) / norm)
  rm(norm)
  
  score = matrix(data=NA, nrow = length(iterUser), ncol=sparseIcm@Dim[2])
  imax = sparseIcm@Dim[2] 
  i2 <- 0
  while (i2 < imax) {  
    i1 = i2
    i2 = i1 + 1000 #taille chunk
    if (i2 > imax) i2 = imax
    cat("\n",i2,"\n")
    sim= as.matrix(t(sparseIcm) %*% sparseIcm[,(i1+1):i2])
    score[,(i1+1):i2] = as.matrix(sparseTrain %*% sim)
    # score[,(i1+1):i2] = sparseTrain %*% sim
  }
  rm(i1, i2, imax, sim)
  
  score = score -1000*as.matrix(sparseTrainR1[iterUser,])
  
  submission=matrix(data = NA, nrow = length(iterUser), ncol = 2)
  colnames(submission) <- c("userId","testItems")
  submission = as.data.frame(submission)
  submission$userId= iterUser
  submission$testItems =do.call(rbind, lapply(1:length(iterUser),FUN=function(i){
    cat(i," ")
    res = list(order(score[i,], decreasing = TRUE)[1:5])
    return(res)
  }))
  return (submission)
}

# Partie Locale
LocalContentBased = function(train){
  train$index = 1:nrow(train)
  usersLocal = train %>% filter(rating>7) %>% group_by(userId) %>% summarise(nb = n()) %>% filter(nb > 1)
  usersLocal = usersLocal$userId
  
  testLocal = train %>% filter(rating>7, userId %in% usersLocal) %>% group_by(userId) %>% sample_frac( size = 0.5)
  trainLocal = train %>% filter(!index %in% testLocal$index)
  
  nbUsers = trainLocal %>% filter(userId %in% usersLocal) %>% group_by(userId) %>% summarise(nb = n())# head(nbUsers)
  
  usersLocal =  sort(usersLocal)
  actual = testLocal %>% group_by(userId) %>% summarise(actual= list(itemId) )
  
  popularity =   trainLocal %>% group_by(itemId) %>% summarise(sum=n())
  popularity = trainLocal %>% filter(userId %in% usersLocal) %>% left_join(popularity) %>% group_by(userId) %>% summarise(popularity = mean(sum))
  
  submission = contentBased(trainLocal, usersLocal)
  submission$actual = actual$actual
  submission$nb = nbUsers$nb
  submission$popularity = popularity$popularity
  submission$MAP10 = do.call(rbind, lapply(1:nrow(submission),FUN=function(i){
    return(mapk(10, submission[i,3], submission[i,2]))
  }))
  return (submission[,4:6])
}

BootstrapContentBased = function(train, n){
  sub = c()
  for ( i in 1:n){
    cat("Etape", i," du bootstrap ContentBased \n")
    sub = rbind(sub,LocalContentBased(train))
  }
  return(sub)
}

LocalURMsvd = function(train, Ndim, nUsers){
  train$index = 1:nrow(train)
  usersLocal = train %>% filter(rating>7) %>% group_by(userId) %>% summarise(nb = n()) %>% filter(nb > 1)
  usersLocal = sample(usersLocal$userId, nUsers)
  
  testLocal = train %>% group_by(itemId) %>% mutate( pop = n()) %>%filter(rating>7, userId %in% usersLocal) %>% group_by(userId) %>% arrange(pop)%>% slice(1:(n()/2))    
  trainLocal = train %>% filter(!index %in% testLocal$index)
  
  nbUsers = trainLocal %>% filter(userId %in% usersLocal) %>% group_by(userId) %>% summarise(nb = n())# head(nbUsers)
  
  usersLocal =  sort(usersLocal)
  actual = testLocal %>% group_by(userId) %>% summarise(actual= list(itemId) )
  
  popularity =   trainLocal %>% group_by(itemId) %>% summarise(sum=n())
  popularity = trainLocal %>% filter(userId %in% usersLocal) %>% left_join(popularity) %>% group_by(userId) %>% summarise(popularity = mean(sum))
  
  submission = URMsvd(trainLocal, usersLocal,950)
  submission$actual = actual$actual
  submission$nb = nbUsers$nb
  submission$popularity = popularity$popularity
  submission$MAP10 = do.call(rbind, lapply(1:nrow(submission),FUN=function(i){
    return(mapk(10, submission[i,3], submission[i,2]))
  }))
  return (submission[,4:6])
}

BootstrapURMsvd = function(train, n, Ndim, nUsers){
  sub = c()
  for ( i in 1:n){
    cat("Etape", i," du bootstrap URMsvd \n")
    sub = rbind(sub,LocalURMsvd(train, Ndim, nUsers))
  }
  return(sub)
}

monPlot = function(submission, name){
  cat("Score : ",mean(submission$MAP10) ,"\n")
  # Model = f(popularity)
  x = unique(trunc(quantile(submission$popularity, probs = seq(from = 0, to=1, length.out=11))))
  submission$cut = as.character(cut(submission$popularity, breaks = x ))
  res = submission %>% group_by(cut) %>% summarise( mean = mean(MAP10))
  res$mean = as.numeric(res$mean) #Otherwise it doesn't work, strange...
  barplot(res$mean, names.arg = as.character(x), xlab = "Mean popularity of items rated", ylab="MAP@5" , main=name)
  
  # Model = f(nb)
  x = unique(trunc(quantile(submission$nb, probs = seq(from = 0, to=1, length.out=21))))
  submission$cut = as.character(cut(submission$nb, breaks = x ))
  res = submission %>% group_by(cut) %>% summarise( mean = mean(MAP10))
  res$mean = as.numeric(res$mean) #Otherwise it doesn't work, strange...
  barplot(res$mean, names.arg = as.character(x), xlab = "Number of items rated", ylab="MAP@5", main=name)
  
  # Model = f(nb*popularity)
  var = submission$nb * submission$popularity
  x = unique(trunc(quantile(var, probs = seq(from = 0, to=1, length.out=11))))
  submission$cut = as.character(cut(var, breaks = x ))
  res = submission %>% group_by(cut) %>% summarise( mean = mean(MAP10))
  res$mean = as.numeric(res$mean) #Otherwise it doesn't work, strange...
  barplot(res$mean, names.arg = as.character(x), main=name)
}
##################################

setwd("C:/Users/Jacques/Documents/polimi/Reccomender")

train=read_csv("train.csv")

iterUser = read_csv("C:/Users/Jacques/Documents/polimi/Reccomender/samplesubmission.csv")
iterUser= iterUser$userId

# Test perso
nUsers = 500
train$index = 1:nrow(train)
usersLocal = train %>% filter(rating>7) %>% group_by(userId) %>% summarise(nb = n()) %>% filter(nb > 1)
usersLocal = sample(usersLocal$userId, nUsers)

testLocal = train %>% group_by(itemId) %>% mutate( pop = n()) %>%filter(rating>7, userId %in% usersLocal) %>% group_by(userId) %>% arrange(pop)%>% slice(1:(n()/2))    




# Chercher la meilleure k SVD
iterDim = seq(from = 200, to=2000, by=200)
test = URMsvdBest3(train, iterDim, 20, 500)
plot(iterDim,test, xlab= "Number of dimensions", ylab = "MAP@5", main="Local SVD's performance")

# Partie locale
sub.BootstrapContentBased = BootstrapContentBased(train, 10)
monPlot(sub.BootstrapContentBased, "Content Based")

sub.BootstrapURMsvd = BootstrapURMsvd(train, 20, 950, 500)
monPlot(sub.BootstrapURMsvd, "SVD")

# Partie online
hybrid = function(train, iterUser, Ndim){
  nbUsers = train %>% filter(userId %in% iterUser) %>% group_by(userId) %>% summarise(nb = n())# head(nbUsers)
  
  popularity = train %>% group_by(itemId) %>% summarise(sum=n())
  popularity = train %>% filter(userId %in% iterUser) %>% left_join(popularity) %>% group_by(userId) %>% summarise(popularity = mean(sum))
  
  submission = data.frame(userId = iterUser)
  submission$nb = nbUsers$nb
  submission$popularity = popularity$popularity
  submission$URMsvd = URMsvd(train, iterUser,Ndim)$testItems
  submission$contentBased = contentBased(train, iterUser)$testItems
  # submission$testItems = submission$contentBased
#   cb = which(submission$popularity < 10 | submission$popularity >29)
#   submission$testItems[cb] = submission$contentBased[cb]
  submission$testItems = do.call(rbind, lapply(1:length(iterUser),FUN=function(i){
    svd = unlist(submission$URMsvd[i])
    cb = unlist(submission$contentBased[i])
    res = c(cb[1:4], svd[1])
    return(paste( res, collapse = " "))
  }))[,]
  return (submission[,c(1,6)])
}

Ndim=950
res = hybrid(train, iterUser, Ndim)
nameFile=paste(c("Hybrid.csv"),collapse="")
write.csv(res,file=nameFile, row.names = FALSE)
