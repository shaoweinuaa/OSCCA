# Lasso regularized Cox regression using leave-one-out (LOO) 
# cross validation (CV).
# During each run of LOO CV, 10-fold CV is performed on training set to 
# select the best model. Then the selected model is applied to the single
# held-out test sample to predict death risk.
rm(list=ls())
library(OIsurv)
library(glmnet)

set.seed(1)
ptm <- proc.time()
#606 194 305 614
mydata = read.table("C:\\Users\\shaowei\\Desktop\\jiaojie\\MICCAI2018\\selectedFeature.txt", header = FALSE,sep=',')
mydata = data.matrix(mydata)
s1 = nrow(mydata)
s2 = ncol(mydata)
mySurv = Surv(mydata[, 1], mydata[, 2]);
x = scale(mydata[, 3:s2])
print(s1)
# leave-one-out CV
group = cbind(numeric(s1))
ind = 1:s1
for(i in 1:s1){
  #cvfit = cv.glmnet(x[ind!=i,], mySurv[ind!=i,], family = "cox", maxit=1000, thresh = 1e-03,alpha=0.5,lambda.min.ratio=0.005)
  #cvfit = cv.glmnet(x[ind!=i,], mySurv[ind!=i,], family = "cox",alpha=0,maxit=5000,nlambda=400,lambda.min.ratio=0.005)
  cvfit = cv.glmnet(x[ind!=i,], mySurv[ind!=i,], family = "cox",maxit=500,lambda.min.ratio=1e-6,nlambda=10)
  preTrain = predict(cvfit, newx = x[ind!=i,], s = 0, type="response")
  print(coef(cvfit,s=0))
  mv = median(preTrain)
  preTest = predict(cvfit, newx = x[ind==i,], s = "lambda.min", type="response")
  if(preTest < mv){
    group[i] = 1
  }else{
    group[i] = 2
  }
  print(i)
}

write.table(group, file = "rank.txt",
            row.names = F, col.names = F, sep="\t")

# logrank
log1 = survdiff(mySurv ~ group,rho = 0)
p = pchisq(log1$chisq, 1, lower.tail=FALSE)
print(p)

# plot KM curve
fit = survfit(mySurv ~ group)
n1 = sum(group==1)
leg1 = paste("Low risk(", n1, ")", sep = "")
n2 = sum(group==2)
leg2 = paste("High risk(", n2, ")", sep = "")
leg3= paste("p=", formatC(p, format="g", digits = 3))

png(filename = "OSCCA.png", width = 5.5, height = 5.5,
    units = "cm", res = 300, pointsize = 7)
plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:3,
     col = 1:3, cex = 0.5)
grid()
legend(x = "bottomright", legend = c(leg1, leg2,leg3), lty = 1:3,
       col = 1:3, cex = 0.65)

dev.off()




