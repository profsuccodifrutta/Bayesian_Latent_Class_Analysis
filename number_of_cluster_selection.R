library(Rmixmod)
# convert all to factors
data$Business = as.factor(data$Business)
data$Shop = as.factor(data$Shop)
data$Upfront = as.factor(data$Upfront)
data$Service = as.factor(data$Service)
data$Customer = as.factor(data$Customer)
data$Quarter = as.factor(data$Quarter)

data = data[,-6]


## select K
nb = 2:30
set.seed(44)
xem <- mixmodCluster(data, nbCluster = nb, 
                     criterion = c("BIC","ICL","NEC"),
                     model=mixmodMultinomialModel(listModels="Binary_pk_Ekj"),
                     seed = 182)

## select K
nbC = NULL
criteria = matrix(NA,length(nb),3)
colnames(criteria) = xem@criterion
for(i in 1:length(nb)){
  criteria[i,] = xem@results[[i]]@criterionValue
  nbC[i] = xem@results[[i]]@nbCluster
}

criteria = criteria[order(nbC),]

# Figure 4
x11()
par(mfrow=c(1,2))
for(j in 1:2){
  plot(nb,criteria[,j],ylab=colnames(criteria)[j],type="b",
       xlab="k",cex.axis=0.8,axes=F)
  axis(1,seq(min(nb),max(nb),1))
  axis(2);  box()
}