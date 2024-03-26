
score = function(y, ychap){
	rmse = sqrt(mean((y - ychap)^2))
	mae = mean(abs(y - ychap))
	EV = (1 - sum((y - ychap)^2) / sum((y - mean(y))^2))*100
	Sc = c(rmse,mae,EV)
	names(Sc) = c("RMSE", "MAE", "EV (%)")
	score = round(Sc,2)
}

###########

scorem = function(y, ypred, noms.model){
  npred = ncol(ypred)
  SC = NULL
  for (j in 1:npred){
	SC = cbind(SC, score(y, ypred[,j]))
  }
  colnames(SC) = noms.model 
  scorem = SC
}
