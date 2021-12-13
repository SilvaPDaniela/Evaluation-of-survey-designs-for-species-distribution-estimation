library("Metrics")

# Read the predicted data in observed locations using the designs obtained from 
# previous R file
# This data needs to contain the coordinates (e.g. longitude and latitude), the 
# observed and predicted values for each species, and the corresponding design.
pred_specie1<-read.csv2("./Github/Data/pred_specie1.csv")
pred_specie2<-read.csv2("./Github/Data/pred_specie2.csv")

# Moreover, the stratum and the corresponding area are needed in the case of 
# comparing designs by stratified mean and stratified variance.

# Create the function to compute the assessment criteria of the selection measures
assessment_function<-function(x,obs_col_name,pred_col_name){
  mae<-mae(x[,obs_col_name],x[,pred_col_name])
  rmse<-rmse(x[,obs_col_name],x[,pred_col_name])
  stratum.table<-data.frame(stratum=names(tapply(x[,pred_col_name],x[,"stratum"],mean)),mean=tapply(x[,pred_col_name],x[,"stratum"],mean),
                       var=tapply(x[,pred_col_name],x[,"stratum"],var),length=tapply(x[,pred_col_name],x[,"stratum"],length))
  stratum.table<-unique(merge(stratum.table,pred_specie1[,c("stratum","Area")],all.x=T))
  str_mean<-sum(stratum.table$Area*stratum.table$mean,na.rm=T)/sum(stratum.table$Area,na.rm=T)
  str_var<-sum(((stratum.table$Area/sum(stratum.table$Area,na.rm=T))^2)*((stratum.table$Area-stratum.table$length)/stratum.table$Area)*(stratum.table$var/stratum.table$length),na.rm=T)
  c(mae,rmse,str_mean,str_var)
}

# Function that allows to apply the previous function to all designs for one specie
assess.each.specie<-function(all.x,obs_col_name,pred_col_name){
  assess.df<-sapply(unique(all.x[,"sample"]),function(s) assessment_function(all.x[which(all.x[,"sample"]==s),],obs_col_name,pred_col_name))
  row.names(assess.df)<-c("MAE","RMSE","STR.MEAN","STR.VAR")
  colnames(assess.df)<-paste("Sample",1:length(unique(all.x[,"sample"])),sep="_")
  assess.df
}

# Apply the previous function for each specie
assess.specie1<-assess.each.specie(pred_specie1,"obs_abund_specie1","pred_abund_specie1")
assess.specie2<-assess.each.specie(pred_specie2,"obs_abund_specie2","pred_abund_specie2")

# Save the results
write.csv2(assess.specie1,file="./Github/assess_specie1.csv")
write.csv2(assess.specie2,file="./Github/assess_specie2.csv")