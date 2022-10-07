
rmse_data <- read.table("../res/tn_RMSE_vs_EOBS_years.csv",sep=";",col.names=c("blabla","years","rmse_era5","rmse_molo","rmse_bola"))
my_min <- min(c(rmse_data$rmse_era5,rmse_data$rmse_molo,rmse_data$rmse_bola))-0.1
my_max <- max(c(rmse_data$rmse_era5,rmse_data$rmse_molo,rmse_data$rmse_bola))+0.1
plot(rmse_data$rmse_era5~rmse_data$years, type="l", col="red",lwd=3,xlab="Years",ylab="RMSE",xaxt="none",ylim=c(my_min,my_max))
points(rmse_data$rmse_era5~rmse_data$years, pch=19, col="red")
lines(rmse_data$rmse_molo~rmse_data$years,col="orange",lwd=3)
points(rmse_data$rmse_molo~rmse_data$years,col="orange",pch=19)
lines(rmse_data$rmse_bola~rmse_data$years,col="blue",lwd=3)
points(rmse_data$rmse_bola~rmse_data$years,col="blue",pch=19)
axis(1, seq(min(rmse_data$years),max(rmse_data$years),2),las=2)
legend("topright",c("ERA5-Land","MOLOCH","BOLAM"),lty=1:1,lwd=c(3,3,3),col=c("red","orange","blue"))

corr_data <- read.table("../res/tn_CORR_vs_EOBS_years.csv",sep=";",col.names=c("blabla","years","corr_era5","corr_molo","corr_bola"))
my_min <- 0.8
my_max <- 1
plot(corr_data$corr_era5~corr_data$years, type="l", col="red",lwd=3,xlab="Years",ylab="Correlation",xaxt="none",ylim=c(my_min,my_max))
lines(corr_data$corr_molo~corr_data$years,col="orange",lwd=3)
lines(corr_data$corr_bola~corr_data$years,col="blue",lwd=3)
axis(1, seq(min(corr_data$years),max(corr_data$years),2),las=2)
legend("topright",c("ERA5-Land","MOLOCH","BOLAM"),lty=1:1,lwd=c(3,3,3),col=c("red","orange","blue"))

bias_data <- read.table("../res/tn_BIAS_vs_EOBS_years.csv",sep=";",col.names=c("blabla","years","bias_era5","bias_molo","bias_bola"))
my_min <- 0
my_max <- 1.5
plot(bias_data$bias_era5~bias_data$years, type="l", col="red",lwd=3,xlab="Years",ylab="Multiplicative bias",xaxt="none",ylim=c(my_min,my_max))
points(bias_data$bias_era5~bias_data$years, pch=19, col="red")
lines(bias_data$bias_molo~bias_data$years,col="orange",lwd=3)
points(bias_data$bias_molo~bias_data$years,col="orange",pch=19)
lines(bias_data$bias_bola~bias_data$years,col="blue",lwd=3)
points(bias_data$bias_bola~bias_data$years,col="blue",pch=19)
axis(1, seq(min(bias_data$years),max(bias_data$years),2),las=2)
legend("topright",c("ERA5-Land","MOLOCH","BOLAM"),lty=1:1,lwd=c(3,3,3),col=c("red","orange","blue"))


