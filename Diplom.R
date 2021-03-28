list.dirs("D:/RData/DIPLOM")
library(reshape)
library(moments)
library(gplots)
library(gclus)
library(gvlma)
library(psych)
library(ggm)
library(sm)
library(car)
library(MASS)
library(fitdistrplus)
library(nortest)

dataFull <- read.csv("D:/RData/DIPLOM/dataFull.csv", header=TRUE, sep=";")
data <- as.data.frame(dataFull[,c("VVP", "SG4Z", "X4BR", "X4BRZ", "SDDN", "SNZP", "DKB", "DMAS", "NALOBR", "INVOSN", "INDPOTR")])
dataPsych <- describe(data)

head(data) #первые 10 сгенерированных значений
str(data)

par(mfrow = c(2,2))



#Формула Стерджесса, для подсчёта числа интервалов
nVVP = 1+3.322*log10(length(data$VVP))
nVVP = ceiling(nVVP)

nSG4Z = 1+3.322*log10(length(data$SG4Z))
nSG4Z = ceiling(nSG4Z)

nX4BR=1+3.322*log10(length(data$X4BR))
nX4BR = ceiling(nX4BR)

nX4BRZ = 1+3.322*log10(length(data$X4BRZ))
nX4BRZ = ceiling(nX4BRZ)

nSDDN = 1+3.322*log10(length(data$SDDN))
nSDDN = ceiling(nSDDN)

nSNZP = 1+3.322*log10(length(data$SNZP))
nSNZP = ceiling(nSNZP)

nDKB = 1+3.322*log10(length(data$DKB))
nDKB = ceiling(nDKB)

nDMAS = 1+3.322*log10(length(data$DMAS))
nDMAS = ceiling(nDMAS)

nNALOBR = 1+3.322*log10(length(data$NALOBR))
nNALOBR = ceiling(nNALOBR)

nINVOSN = 1+3.322*log10(length(data$INVOSN))
nINVOSN = ceiling(nINVOSN)

nINDPOTR = 1+3.322*log10(length(data$INDPOTR))
nINDPOTR = ceiling(nINDPOTR)


#Гистограммы. Частота
hist(data$VVP, freq=TRUE, col="blue", xlab="VVP", main="Гистограмма. Частота")
hist(data$SG4Z, freq=TRUE, col="blue", xlab="SG4Z", main="Гистограмма. Частота")
hist(data$X4BR, freq=TRUE, col="blue", xlab="X4BR", main="Гистограмма. Частота")
hist(data$X4BRZ, freq=TRUE, col="blue", xlab="X4BRZ", main="Гистограмма. Частота")
hist(data$SDDN, freq=TRUE, col="blue", xlab="SDDN", main="Гистограмма. Частота")
hist(data$SNZP, freq=TRUE, col="blue", xlab="SNZP", main="Гистограмма. Частота")
hist(data$DKB, freq=TRUE, col="blue", xlab="DKB", main="Гистограмма. Частота")
hist(data$DMAS, freq=TRUE, col="blue", xlab="DMAS", main="Гистограмма. Частота")
hist(data$NALOBR, freq=TRUE, col="blue", xlab="NALOBR", main="Гистограмма. Частота")
hist(data$INVOSN, freq=TRUE, col="blue", xlab="INVOSN", main="Гистограмма. Частота")
hist(data$INDPOTR, freq=TRUE, col="blue", xlab="INDPOTR", main="Гистограмма. Частота")

#Гистограммы и Ядерная плотность
#dev.new()
#VVP
h = hist(
  data$VVP,
  breaks = nVVP,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "VVP"
)
xfit=seq(min(data$VVP), max(data$VVP), length=15)
yfit=dnorm(xfit, mean=mean(data$VVP), sd=sd(data$VVP))
yfit=yfit*diff(h$mids[1:2])*length(data$VVP)
lines(xfit, yfit, col="red", lwd=2)

#SG4Z
h = hist(
  data$SG4Z,
  breaks = nSG4Z,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "SG4Z"
)
xfit=seq(min(data$SG4Z), max(data$SG4Z), length=15)
yfit=dnorm(xfit, mean=mean(data$SG4Z), sd=sd(data$SG4Z))
yfit=yfit*diff(h$mids[1:2])*length(data$SG4Z)
lines(xfit, yfit, col="red", lwd=2)

#X4BR
h = hist(
  data$X4BR,
  breaks = nX4BR,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "X4BR"
)
xfit=seq(min(data$X4BR), max(data$X4BR), length=15)
yfit=dnorm(xfit, mean=mean(data$X4BR), sd=sd(data$X4BR))
yfit=yfit*diff(h$mids[1:2])*length(data$X4BR)
lines(xfit, yfit, col="red", lwd=2)

#X4BRZ
h = hist(
  data$X4BRZ,
  breaks = nX4BRZ,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "X4BRZ"
)
xfit=seq(min(data$X4BRZ), max(data$X4BRZ), length=15)
yfit=dnorm(xfit, mean=mean(data$X4BRZ), sd=sd(data$X4BRZ))
yfit=yfit*diff(h$mids[1:2])*length(data$X4BRZ)
lines(xfit, yfit, col="red", lwd=2)

#SDDN
h = hist(
  data$SDDN,
  breaks = nSDDN,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "SDDN"
)
xfit=seq(min(data$SDDN), max(data$SDDN), length=15)
yfit=dnorm(xfit, mean=mean(data$SDDN), sd=sd(data$SDDN))
yfit=yfit*diff(h$mids[1:2])*length(data$SDDN)
lines(xfit, yfit, col="red", lwd=2)

#SNZP
h = hist(
  data$SNZP,
  breaks = nSNZP,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "SNZP"
)
xfit=seq(min(data$SNZP), max(data$SNZP), length=15)
yfit=dnorm(xfit, mean=mean(data$SNZP), sd=sd(data$SNZP))
yfit=yfit*diff(h$mids[1:2])*length(data$SNZP)
lines(xfit, yfit, col="red", lwd=2)

#DKB
h = hist(
  data$DKB,
  breaks = nDKB,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "DKB"
)
xfit=seq(min(data$DKB), max(data$DKB), length=15)
yfit=dnorm(xfit, mean=mean(data$DKB), sd=sd(data$DKB))
yfit=yfit*diff(h$mids[1:2])*length(data$DKB)
lines(xfit, yfit, col="red", lwd=2)

#DMAS
h = hist(
  data$DMAS,
  breaks = nDMAS,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "DMAS"
)
xfit=seq(min(data$DMAS), max(data$DMAS), length=15)
yfit=dnorm(xfit, mean=mean(data$DMAS), sd=sd(data$DMAS))
yfit=yfit*diff(h$mids[1:2])*length(data$DMAS)
lines(xfit, yfit, col="red", lwd=2)

#NALOBR
h = hist(
  data$NALOBR,
  breaks = nNALOBR,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "NALOBR"
)
xfit=seq(min(data$NALOBR), max(data$NALOBR), length=15)
yfit=dnorm(xfit, mean=mean(data$NALOBR), sd=sd(data$NALOBR))
yfit=yfit*diff(h$mids[1:2])*length(data$NALOBR)
lines(xfit, yfit, col="red", lwd=2)

#INVOSN
h = hist(
  data$INVOSN,
  breaks = nINVOSN,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "INVOSN"
)
xfit=seq(min(data$INVOSN), max(data$INVOSN), length=15)
yfit=dnorm(xfit, mean=mean(data$INVOSN), sd=sd(data$INVOSN))
yfit=yfit*diff(h$mids[1:2])*length(data$INVOSN)
lines(xfit, yfit, col="red", lwd=2)

#INDPOTR
h = hist(
  data$INDPOTR,
  breaks = nINDPOTR,
  main = "Гистограмма с числом интервалов n и с кривой нормального распределения",
  col = "blue",
  xlab = "INDPOTR"
)
xfit=seq(min(data$INDPOTR), max(data$INDPOTR), length=15)
yfit=dnorm(xfit, mean=mean(data$INDPOTR), sd=sd(data$INDPOTR))
yfit=yfit*diff(h$mids[1:2])*length(data$INDPOTR)
lines(xfit, yfit, col="red", lwd=2)


#Визуализация корреляций
scatterplotMatrix(data, spread=FALSE, lty.smooth=2, main="Диаграмма рассеяния")

#Проверка статистической значимости корреляций
cor.test(data$VVP, data$SG4Z, alternative="two.side", method="pearson")
cor.test(data$VVP, data$X4BR, alternative="two.side", method="pearson")
cor.test(data$VVP, data$X4BRZ, alternative="two.side", method="pearson")
cor.test(data$VVP, data$SDDN, alternative="two.side", method="pearson")
cor.test(data$VVP, data$SNZP, alternative="two.side", method="pearson")
cor.test(data$VVP, data$DKB, alternative="two.side", method="pearson")
cor.test(data$VVP, data$DMAS, alternative="two.side", method="pearson")
cor.test(data$VVP, data$NALOBR, alternative="two.side", method="pearson")
cor.test(data$VVP, data$INVOSN, alternative="two.side", method="pearson")
cor.test(data$VVP, data$INDPOTR, alternative="two.side", method="pearson")

#Оценка параметров распределений
graph_distr <- function (x, pc, pd, main_name="")
{
  op <- par(mfrow=c(2,2), pty="s")
  par(mfrow=c(2,2))
  mn <- paste(c("Эмпир-я КФР и ", main_name))
  plot(x,pc, type="l",col="red", lwd=2, main=mn)
  plot(ecdf(x),add=TRUE)
  mn <- paste(c("Эмпир-я ФПР и ", main_name))
  plot(density(x), lwd = 2, col="blue", main=mn)
  lines(x, pd, col="red",lwd = 2)
  par(op)
}

par(mfrow=c(2,2))

##Оценка параметров нормального распределения
(dof_VVP = fitdistr(data$VVP,"normal"))
ep1_VVP=dof_VVP$estimate[1];
ep2_VVP=dof_VVP$estimate[2];
graph_distr(data$VVP, pnorm(data$VVP,mean=ep1_VVP,sd=ep2_VVP), dnorm(data$VVP,mean=ep1_VVP,sd=ep2_VVP), "норм. р., VVP")

(dof_SG4Z = fitdistr(data$SG4Z,"normal"))
ep1_SG4Z=dof_SG4Z$estimate[1];
ep2_SG4Z=dof_SG4Z$estimate[2];
graph_distr(data$SG4Z, pnorm(data$SG4Z,mean=ep1_SG4Z,sd=ep2_SG4Z), dnorm(data$SG4Z,mean=ep1_SG4Z,sd=ep2_SG4Z), "норм. р., SG4Z")

(dof_X4BR = fitdistr(data$X4BR,"normal"))
ep1_X4BR=dof_X4BR$estimate[1];
ep2_X4BR=dof_X4BR$estimate[2];
graph_distr(data$X4BR, pnorm(data$X4BR,mean=ep1_X4BR,sd=ep2_X4BR), dnorm(data$X4BR,mean=ep1_X4BR,sd=ep2_X4BR), "норм. р., X4BR")

(dof_X4BRZ = fitdistr(data$X4BRZ,"normal"))
ep1_X4BRZ=dof_X4BRZ$estimate[1];
ep2_X4BRZ=dof_X4BRZ$estimate[2];
graph_distr(data$X4BRZ, pnorm(data$X4BRZ,mean=ep1_X4BRZ,sd=ep2_X4BRZ), dnorm(data$X4BRZ,mean=ep1_X4BRZ,sd=ep2_X4BRZ), "норм. р., X4BRZ")

(dof_SDDN = fitdistr(data$SDDN,"normal"))
ep1_SDDN=dof_SDDN$estimate[1];
ep2_SDDN=dof_SDDN$estimate[2];
graph_distr(data$SDDN, pnorm(data$SDDN,mean=ep1_SDDN,sd=ep2_SDDN), dnorm(data$SDDN,mean=ep1_SDDN,sd=ep2_SDDN), "норм. р., SDDN")

(dof_SNZP = fitdistr(data$SNZP,"normal"))
ep1_SNZP=dof_SNZP$estimate[1];
ep2_SNZP=dof_SNZP$estimate[2];
graph_distr(data$SNZP, pnorm(data$SNZP,mean=ep1_SNZP,sd=ep2_SNZP), dnorm(data$SNZP,mean=ep1_SNZP,sd=ep2_SNZP), "норм. р., SNZP")

(dof_DKB = fitdistr(data$DKB,"normal"))
ep1_DKB=dof_DKB$estimate[1];
ep2_DKB=dof_DKB$estimate[2];
graph_distr(data$DKB, pnorm(data$DKB,mean=ep1_DKB,sd=ep2_DKB), dnorm(data$DKB,mean=ep1_DKB,sd=ep2_DKB), "норм. р., DKB")

#DMAS
(dof_DMAS = fitdistr(data$DMAS,"normal"))
ep1_DMAS=dof_DMAS$estimate[1];
ep2_DMAS=dof_DMAS$estimate[2];
graph_distr(data$DMAS, pnorm(data$DMAS,mean=ep1_DMAS,sd=ep2_DMAS), dnorm(data$DMAS,mean=ep1_DMAS,sd=ep2_DMAS), "норм. р., DMAS")

#NALOBR
(dof_NALOBR = fitdistr(data$NALOBR,"normal"))
ep1_NALOBR=dof_NALOBR$estimate[1];
ep2_NALOBR=dof_NALOBR$estimate[2];
graph_distr(data$NALOBR, pnorm(data$NALOBR,mean=ep1_NALOBR,sd=ep2_NALOBR), dnorm(data$NALOBR,mean=ep1_NALOBR,sd=ep2_NALOBR), "норм. р., NALOBR")

#INVOSN
(dof_INVOSN = fitdistr(data$INVOSN,"normal"))
ep1_INVOSN=dof_INVOSN$estimate[1];
ep2_INVOSN=dof_INVOSN$estimate[2];
graph_distr(data$INVOSN, pnorm(data$INVOSN,mean=ep1_INVOSN,sd=ep2_INVOSN), dnorm(data$INVOSN,mean=ep1_INVOSN,sd=ep2_INVOSN), "норм. р., INVOSN")

#INDPOTR
(dof_INDPOTR = fitdistr(data$INDPOTR,"normal"))
ep1_INDPOTR=dof_INDPOTR$estimate[1];
ep2_INDPOTR=dof_INDPOTR$estimate[2];
graph_distr(data$INDPOTR, pnorm(data$INDPOTR,mean=ep1_INDPOTR,sd=ep2_INDPOTR), dnorm(data$INDPOTR,mean=ep1_INDPOTR,sd=ep2_INDPOTR), "норм. р., INDPOTR")


##Оценка параметров лог-нормального распределения
(dof_VVP_lg = fitdistr(data$VVP,"log-normal"))
ep1_VVP_lg=dof_VVP_lg$estimate[1];
ep2_VVP_lg=dof_VVP_lg$estimate[2];
graph_distr(data$VVP, plnorm(data$VVP,meanlog=ep1_VVP_lg,sdlog=ep2_VVP_lg), dlnorm(data$VVP,meanlog=ep1_VVP_lg,sdlog=ep2_VVP_lg), "лог-норм. р., VVP")

(dof_SG4Z_lg = fitdistr(data$SG4Z,"log-normal"))
ep1_SG4Z_lg=dof_SG4Z_lg$estimate[1];
ep2_SG4Z_lg=dof_SG4Z_lg$estimate[2];
graph_distr(data$SG4Z, plnorm(data$SG4Z,meanlog=ep1_SG4Z_lg,sdlog=ep2_SG4Z_lg), dlnorm(data$SG4Z,meanlog=ep1_SG4Z_lg,sdlog=ep2_SG4Z_lg), "лог-норм. р., SG4Z")

(dof_X4BR_lg = fitdistr(data$X4BR,"log-normal"))
ep1_X4BR_lg=dof_X4BR_lg$estimate[1];
ep2_X4BR_lg=dof_X4BR_lg$estimate[2];
graph_distr(data$X4BR, plnorm(data$X4BR,meanlog=ep1_X4BR_lg,sdlog=ep2_X4BR_lg), dlnorm(data$X4BR,meanlog=ep1_X4BR_lg,sdlog=ep2_X4BR_lg), "лог-норм. р., X4BR")

(dof_X4BRZ_lg = fitdistr(data$X4BRZ,"log-normal"))
ep1_X4BRZ_lg=dof_X4BRZ_lg$estimate[1];
ep2_X4BRZ_lg=dof_X4BRZ_lg$estimate[2];
graph_distr(data$X4BRZ, plnorm(data$X4BRZ,meanlog=ep1_X4BRZ_lg,sdlog=ep2_X4BRZ_lg), dlnorm(data$X4BRZ,meanlog=ep1_X4BRZ_lg,sdlog=ep2_X4BRZ_lg), "лог-норм. р., X4BRZ")

(dof_SDDN_lg = fitdistr(data$SDDN,"log-normal"))
ep1_SDDN_lg=dof_SDDN_lg$estimate[1];
ep2_SDDN_lg=dof_SDDN_lg$estimate[2];
graph_distr(data$SDDN, plnorm(data$SDDN,meanlog=ep1_SDDN_lg,sdlog=ep2_SDDN_lg), dlnorm(data$SDDN,meanlog=ep1_SDDN_lg,sdlog=ep2_SDDN_lg), "лог-норм. р., SDDN")

(dof_SNZP_lg = fitdistr(data$SNZP,"log-normal"))
ep1_SNZP_lg=dof_SNZP_lg$estimate[1];
ep2_SNZP_lg=dof_SNZP_lg$estimate[2];
graph_distr(data$SNZP, plnorm(data$SNZP,meanlog=ep1_SNZP_lg,sdlog=ep2_SNZP_lg), dlnorm(data$SNZP,meanlog=ep1_SNZP_lg,sdlog=ep2_SNZP_lg), "лог-норм. р., SNZP")

(dof_DKB_lg = fitdistr(data$DKB,"log-normal"))
ep1_DKB_lg=dof_DKB_lg$estimate[1];
ep2_DKB_lg=dof_DKB_lg$estimate[2];
graph_distr(data$DKB, plnorm(data$DKB,meanlog=ep1_DKB_lg,sdlog=ep2_DKB_lg), dlnorm(data$DKB,meanlog=ep1_DKB_lg,sdlog=ep2_DKB_lg), "лог-норм. р., DKB")

(dof_DMAS_lg = fitdistr(data$DMAS,"log-normal"))
ep1_DMAS_lg=dof_DMAS_lg$estimate[1];
ep2_DMAS_lg=dof_DMAS_lg$estimate[2];
graph_distr(data$DMAS, plnorm(data$DMAS,meanlog=ep1_DMAS_lg,sdlog=ep2_DMAS_lg), dlnorm(data$DMAS,meanlog=ep1_DMAS_lg,sdlog=ep2_DMAS_lg), "лог-норм. р., DMAS")

(dof_NALOBR_lg = fitdistr(data$NALOBR,"log-normal"))
ep1_NALOBR_lg=dof_NALOBR_lg$estimate[1];
ep2_NALOBR_lg=dof_NALOBR_lg$estimate[2];
graph_distr(data$NALOBR, plnorm(data$NALOBR,meanlog=ep1_NALOBR_lg,sdlog=ep2_NALOBR_lg), dlnorm(data$NALOBR,meanlog=ep1_NALOBR_lg,sdlog=ep2_NALOBR_lg), "лог-норм. р., NALOBR")

(dof_INVOSN_lg = fitdistr(data$INVOSN,"log-normal"))
ep1_INVOSN_lg=dof_INVOSN_lg$estimate[1];
ep2_INVOSN_lg=dof_INVOSN_lg$estimate[2];
graph_distr(data$INVOSN, plnorm(data$INVOSN,meanlog=ep1_INVOSN_lg,sdlog=ep2_INVOSN_lg), dlnorm(data$INVOSN,meanlog=ep1_INVOSN_lg,sdlog=ep2_INVOSN_lg), "лог-норм. р., INVOSN")

(dof_INDPOTR_lg = fitdistr(data$INDPOTR,"log-normal"))
ep1_INDPOTR_lg=dof_INDPOTR_lg$estimate[1];
ep2_INDPOTR_lg=dof_INDPOTR_lg$estimate[2];
graph_distr(data$INDPOTR, plnorm(data$INDPOTR,meanlog=ep1_INDPOTR_lg,sdlog=ep2_INDPOTR_lg), dlnorm(data$INDPOTR,meanlog=ep1_INDPOTR_lg,sdlog=ep2_INDPOTR_lg), "лог-норм. р., INDPOTR")


##Оценка параметров распределения Вейбулла
(dof_VVP_w = fitdist(data$VVP,"weibull"))
ep1_VVP_w=dof_VVP_w$estimate[1];
ep2_VVP_w=dof_VVP_w$estimate[2];
graph_distr(data$VVP, pweibull(data$VVP,scale=ep1_VVP_w,shape=ep2_VVP_w), dweibull(data$VVP,scale=ep1_VVP_w,shape=ep2_VVP_w),"р.Вейбулла, VVP")

(dof_SG4Z_w = fitdist(data$SG4Z,"weibull"))
ep1_SG4Z_w=dof_SG4Z_w$estimate[1];
ep2_SG4Z_w=dof_SG4Z_w$estimate[2];
graph_distr(data$SG4Z, pweibull(data$SG4Z,scale=ep1_SG4Z_w,shape=ep2_SG4Z_w), dweibull(data$SG4Z,scale=ep1_SG4Z_w,shape=ep2_SG4Z_w), "р.Вейбулла, SG4Z")

(dof_X4BR_w = fitdist(data$X4BR,"weibull"))
ep1_X4BR_w=dof_X4BR_w$estimate[1];
ep2_X4BR_w=dof_X4BR_w$estimate[2];
graph_distr(data$X4BR, pweibull(data$X4BR,scale=ep1_X4BR_w,shape=ep2_X4BR_w), dweibull(data$X4BR,scale=ep1_X4BR_w,shape=ep2_X4BR_w), "р.Вейбулла, X4BR")

(dof_X4BRZ_w = fitdist(data$X4BRZ,"weibull"))
ep1_X4BRZ_w=dof_X4BRZ_w$estimate[1];
ep2_X4BRZ_w=dof_X4BRZ_w$estimate[2];
graph_distr(data$X4BRZ, pweibull(data$X4BRZ,scale=ep1_X4BRZ_w,shape=ep2_X4BRZ_w), dweibull(data$X4BRZ,scale=ep1_X4BRZ_w,shape=ep2_X4BRZ_w), "р.Вейбулла, X4BRZ")

(dof_SDDN_w = fitdist(data$SDDN,"weibull"))
ep1_SDDN_w=dof_SDDN_w$estimate[1];
ep2_SDDN_w=dof_SDDN_w$estimate[2];
graph_distr(data$SDDN, pweibull(data$SDDN,scale=ep1_SDDN_w,shape=ep2_SDDN_w), dweibull(data$SDDN,scale=ep1_SDDN_w,shape=ep2_SDDN_w), "р.Вейбулла, SDDN")

(dof_SNZP_w = fitdist(data$SNZP,"weibull"))
ep1_SNZP_w=dof_SNZP_w$estimate[1];
ep2_SNZP_w=dof_SNZP_w$estimate[2];
graph_distr(data$SNZP, pweibull(data$SNZP,scale=ep1_SNZP_w,shape=ep2_SNZP_w), dweibull(data$SNZP,scale=ep1_SNZP_w,shape=ep2_SNZP_w), "р.Вейбулла, SNZP")

(dof_DKB_w = fitdist(data$DKB,"weibull"))
ep1_DKB_w=dof_DKB_w$estimate[1];
ep2_DKB_w=dof_DKB_w$estimate[2];
graph_distr(data$DKB, pweibull(data$DKB,scale=ep1_DKB_w,shape=ep2_DKB_w), dweibull(data$DKB,scale=ep1_DKB_w,shape=ep2_DKB_w), "р.Вейбулла, DKB")

(dof_DMAS_w = fitdist(data$DMAS,"weibull"))
ep1_DMAS_w=dof_DMAS_w$estimate[1];
ep2_DMAS_w=dof_DMAS_w$estimate[2];
graph_distr(data$DMAS, pweibull(data$DMAS,scale=ep1_DMAS_w,shape=ep2_DMAS_w), dweibull(data$DMAS,scale=ep1_DMAS_w,shape=ep2_DMAS_w), "р.Вейбулла, DMAS")

(dof_NALOBR_w = fitdist(data$NALOBR,"weibull"))
ep1_NALOBR_w=dof_NALOBR_w$estimate[1];
ep2_NALOBR_w=dof_NALOBR_w$estimate[2];
graph_distr(data$NALOBR, pweibull(data$NALOBR,scale=ep1_NALOBR_w,shape=ep2_NALOBR_w), dweibull(data$NALOBR,scale=ep1_NALOBR_w,shape=ep2_NALOBR_w), "р.Вейбулла, NALOBR")

(dof_INVOSN_w = fitdist(data$INVOSN,"weibull"))
ep1_INVOSN_w=dof_INVOSN_w$estimate[1];
ep2_INVOSN_w=dof_INVOSN_w$estimate[2];
graph_distr(data$INVOSN, pweibull(data$INVOSN,scale=ep1_INVOSN_w,shape=ep2_INVOSN_w), dweibull(data$INVOSN,scale=ep1_INVOSN_w,shape=ep2_INVOSN_w), "р.Вейбулла, INVOSN")

(dof_INDPOTR_w = fitdist(data$INDPOTR,"weibull"))
ep1_INDPOTR_w=dof_INDPOTR_w$estimate[1];
ep2_INDPOTR_w=dof_INDPOTR_w$estimate[2];
graph_distr(data$INDPOTR, pweibull(data$INDPOTR,scale=ep1_INDPOTR_w,shape=ep2_INDPOTR_w), dweibull(data$INDPOTR,scale=ep1_INDPOTR_w,shape=ep2_INDPOTR_w), "р.Вейбулла, INDPOTR")


#Квантиль-квантильные графики
x_VVP <- data$VVP;
cvm.test(x_VVP);
qqnorm(x_VVP); qqline(x_VVP);
z_VVP <- (x_VVP - mean(x_VVP))/sqrt(var(x_VVP)) # Стандартизация выборки
x_VVP.qq <- qqnorm(z_VVP, plot.it = FALSE)
x_VVP.qq <- lapply(x_VVP.qq, sort)
plot(x_VVP.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, VVP")

x_SG4Z <- data$SG4Z;
qqnorm(x_SG4Z); qqline(x_SG4Z);
z_SG4Z <- (x_SG4Z - mean(x_SG4Z))/sqrt(var(x_SG4Z)) # Стандартизация выборки
x_SG4Z.qq <- qqnorm(z_SG4Z, plot.it = FALSE)
x_SG4Z.qq <- lapply(x_SG4Z.qq, sort)
plot(x_SG4Z.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, SG4Z")

x_X4BR <- data$X4BR;
cvm.test(x_X4BR);
qqnorm(x_X4BR); qqline(x_X4BR);
z_X4BR <- (x_X4BR - mean(x_X4BR))/sqrt(var(x_X4BR)) # Стандартизация выборки
x_X4BR.qq <- qqnorm(z_X4BR, plot.it = FALSE)
x_X4BR.qq <- lapply(x_X4BR.qq, sort)
plot(x_X4BR.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, X4BR")

x_X4BRZ <- data$X4BRZ;
qqnorm(x_X4BRZ); qqline(x_X4BRZ);
z_X4BRZ <- (x_X4BRZ - mean(x_X4BRZ))/sqrt(var(x_X4BRZ)) # Стандартизация выборки
x_X4BRZ.qq <- qqnorm(z_X4BRZ, plot.it = FALSE)
x_X4BRZ.qq <- lapply(x_X4BRZ.qq, sort)
plot(x_X4BRZ.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, X4BRZ")

x_SDDN <- data$SDDN;
cvm.test(x_SDDN);
qqnorm(x_SDDN); qqline(x_SDDN);
z_SDDN <- (x_SDDN - mean(x_SDDN))/sqrt(var(x_SDDN)) # Стандартизация выборки
x_SDDN.qq <- qqnorm(z_SDDN, plot.it = FALSE)
x_SDDN.qq <- lapply(x_SDDN.qq, sort)
plot(x_SDDN.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, SDDN")

x_SNZP <- data$SNZP;
qqnorm(x_SNZP); qqline(x_SNZP);
z_SNZP <- (x_SNZP - mean(x_SNZP))/sqrt(var(x_SNZP)) # Стандартизация выборки
x_SNZP.qq <- qqnorm(z_SNZP, plot.it = FALSE)
x_SNZP.qq <- lapply(x_SNZP.qq, sort)
plot(x_SNZP.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, SNZP")

x_DKB <- data$DKB;
qqnorm(x_DKB); qqline(x_DKB);
z_DKB <- (x_DKB - mean(x_DKB))/sqrt(var(x_DKB)) # Стандартизация выборки
x_DKB.qq <- qqnorm(z_DKB, plot.it = FALSE)
x_DKB.qq <- lapply(x_DKB.qq, sort)
plot(x_DKB.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, DKB")

x_DMAS <- data$DMAS;
qqnorm(x_DMAS); qqline(x_DMAS);
z_DMAS <- (x_DMAS - mean(x_DMAS))/sqrt(var(x_DMAS)) # Стандартизация выборки
x_DMAS.qq <- qqnorm(z_DMAS, plot.it = FALSE)
x_DMAS.qq <- lapply(x_DMAS.qq, sort)
plot(x_DMAS.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, DMAS")

x_NALOBR <- data$NALOBR;
qqnorm(x_NALOBR); qqline(x_NALOBR);
z_NALOBR <- (x_NALOBR - mean(x_NALOBR))/sqrt(var(x_NALOBR)) # Стандартизация выборки
x_NALOBR.qq <- qqnorm(z_NALOBR, plot.it = FALSE)
x_NALOBR.qq <- lapply(x_NALOBR.qq, sort)
plot(x_NALOBR.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, NALOBR")

x_INVOSN <- data$INVOSN;
qqnorm(x_INVOSN); qqline(x_INVOSN);
z_INVOSN <- (x_INVOSN - mean(x_INVOSN))/sqrt(var(x_INVOSN)) # Стандартизация выборки
x_INVOSN.qq <- qqnorm(z_INVOSN, plot.it = FALSE)
x_INVOSN.qq <- lapply(x_INVOSN.qq, sort)
plot(x_INVOSN.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, INVOSN")

x_INDPOTR <- data$INDPOTR;
qqnorm(x_INDPOTR); qqline(x_INDPOTR);
z_INDPOTR <- (x_INDPOTR - mean(x_INDPOTR))/sqrt(var(x_INDPOTR)) # Стандартизация выборки
x_INDPOTR.qq <- qqnorm(z_INDPOTR, plot.it = FALSE)
x_INDPOTR.qq <- lapply(x_INDPOTR.qq, sort)
plot(x_INDPOTR.qq, ylim = c(-2, 5), ylab = "Z-статистики выборки", xlab = "Квантили НР, INDPOTR")


#Функции плотности распределения и Тесты на нормальность Шапиро-Уилкса
sm.density(data$VVP, model = "Normal", xlab="VVP, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$VVP)

sm.density(data$SG4Z, model = "Normal", xlab="SG4Z, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$SG4Z)

sm.density(data$X4BR, model = "Normal", xlab="X4BR, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$X4BR)

sm.density(data$X4BRZ, model = "Normal", xlab="X4BRZ, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$X4BRZ)

sm.density(data$SDDN, model = "Normal", xlab="SDDN, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$SDDN)

sm.density(data$SNZP, model = "Normal", xlab="SNZP, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$SNZP)

sm.density(data$DKB, model = "Normal", xlab="DKB, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$DKB)

sm.density(data$DMAS, model = "Normal", xlab="DMAS, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$DMAS)

sm.density(data$NALOBR, model = "Normal", xlab="NALOBR, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$NALOBR)

sm.density(data$INVOSN, model = "Normal", xlab="INVOSN, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$INVOSN)

sm.density(data$INDPOTR, model = "Normal", xlab="INDPOTR, Имитированная выборка", ylab = "Функция плотности распределения")
shapiro.test(data$INDPOTR)




names(dataPsych)[1] <- "ID"
names(dataPsych)[2] <- "Кол-во"
names(dataPsych)[3] <- "Среднее"
names(dataPsych)[4] <- "Станд. Отклонение"
names(dataPsych)[5] <- "Медиана"
names(dataPsych)[6] <- "Отсечение"
names(dataPsych)[7] <- "Среднее станд. отклонение"
names(dataPsych)[8] <- "Мин"
names(dataPsych)[9] <- "Макс"
names(dataPsych)[10] <- "Разброс"
names(dataPsych)[11] <- "Скос"
names(dataPsych)[12] <- "Эксцесс"
names(dataPsych)[13] <- "Средняя квадратическая ошибка"




#######Модель №1#######

statesOP1 <- as.data.frame(dataFull[,c("VVP","SG4Z","X4BR")])

head(statesOP1)

str(statesOP1)

#Регрессия OP1

OP1 <- lm(VVP ~ SG4Z + X4BR, data = statesOP1)

#Кореляция + Матрица диаграмм рассеяния OP1

cor(statesOP1)

scatterplotMatrix(statesOP1, spread = FALSE, lty.smooth = 2,main = "Матрица диаграмм рассеяния OP1")

#Вывод данных OP1

summary(OP1)

confint(OP1)

#Стандартный подход OP1

par(mfrow = c(2,2))

plot(OP1)

#Нормальность OP1

qqPlot(OP1, labels = row.names(statesOP1),id.method = "identify",simulate = TRUE, main = "Q-Q Plot OP1")

residplot <- function(OP1, nbreaks=10){
  
  z <- rstudent(OP1)
  
  hist(z, breaks = nbreaks, freq = FALSE, xlab = "Остатки Стьюдента", main = "Распределение остатков OP1")
  
  rug(jitter(z), col = "brown")
  
  curve(dnorm(x, mean = mean(z), sd = sd(z)), add = TRUE, col = "blue", lwd = 2)
  
  lines(density(z)$x, density(z)$y, col = "red", lwd = 2, lty = 2)
  
  legend("topright",legend = c( "Кривая нормального распределения",
                                
                                "Ядерная оценка плотности"),
         
         lty = 1:2, col = c("blue","red"), cex = .7)
  
}

residplot(OP1)

#Независимость остатков OP1

durbinWatsonTest(OP1)

#Линейность OP1

crPlots(OP1)

#Гомоскедастичность OP1

ncvTest(OP1)

spreadLevelPlot(OP1)

#Общая проверка выполнения требований OP1

gvmodelOP1 <- gvlma(OP1)

summary(gvmodelOP1)

#Мультиколлинеарность OP1

vif(OP1)

sqrt(vif(OP1)) > 2 # проблема?

#######Модель №2#######

statesOP2 <- as.data.frame(dataFull[,c("VVP","X4BRZ","SNZP")])

head(statesOP2)

str(statesOP2)

#Регрессия OP2

OP2 <- lm(VVP ~ X4BRZ + SNZP, data = statesOP2)

#Кореляция + Матрица диаграмм рассеяния OP2

cor(statesOP2)

scatterplotMatrix(statesOP2, spread = FALSE, lty.smooth = 2,main = "Матрица диаграмм рассеяния OP2")

#Вывод данных OP2

summary(OP2)

confint(OP2)

#Стандартный подход OP2

par(mfrow = c(2,2))

plot(OP2)

#Нормальность OP2

qqPlot(OP2, labels = row.names(statesOP2),id.method = "identify",simulate = TRUE, main = "Q-Q Plot OP2")

residplot <- function(OP2, nbreaks=10){
  
  z <- rstudent(OP2)
  
  hist(z, breaks = nbreaks, freq = FALSE, xlab = "Остатки Стьюдента", main = "Распределение остатков OP2")
  
  rug(jitter(z), col = "brown")
  
  curve(dnorm(x, mean = mean(z), sd = sd(z)), add = TRUE, col = "blue", lwd = 2)
  
  lines(density(z)$x, density(z)$y, col = "red", lwd = 2, lty = 2)
  
  legend("topright",legend = c( "Кривая нормального распределения",
                                
                                "Ядерная оценка плотности"),
         
         lty = 1:2, col = c("blue","red"), cex = .7)
  
}

residplot(OP2)

#Независимость остатков OP2

durbinWatsonTest(OP2)

#Линейность OP2

crPlots(OP2)

#Гомоскедастичность OP2

ncvTest(OP2)

spreadLevelPlot(OP2)

#Общая проверка выполнения требований OP2

gvmodelOP2 <- gvlma(OP2)

summary(gvmodelOP2)

#Мультиколлинеарность OP2

vif(OP2)

sqrt(vif(OP2)) > 2 # проблема?

#######Модель №3#######

statesOP3 <- as.data.frame(dataFull[,c("VVP","X4BRZ","SNZP","DKB")])

head(statesOP3)

str(statesOP3)

#Регрессия OP3

OP3 <- lm(VVP ~ X4BRZ + SNZP + DKB, data = statesOP3)

#Кореляция + Матрица диаграмм рассеяния OP3

cor(statesOP3)

scatterplotMatrix(statesOP3, spread = FALSE, lty.smooth = 2,main = "Матрица диаграмм рассеяния OP3")

#Вывод данных OP3

summary(OP3)

confint(OP3)

#Стандартный подход OP3

par(mfrow = c(2,2))

plot(OP3)

#Нормальность OP3

qqPlot(OP3, labels = row.names(statesOP3),id.method = "identify",simulate = TRUE, main = "Q-Q Plot OP3")

residplot <- function(OP3, nbreaks=10){
  
  z <- rstudent(OP3)
  
  hist(z, breaks = nbreaks, freq = FALSE, xlab = "Остатки Стьюдента", main = "Распределение остатков OP3")
  
  rug(jitter(z), col = "brown")
  
  curve(dnorm(x, mean = mean(z), sd = sd(z)), add = TRUE, col = "blue", lwd = 2)
  
  lines(density(z)$x, density(z)$y, col = "red", lwd = 2, lty = 2)
  
  legend("topright",legend = c( "Кривая нормального распределения",
                                
                                "Ядерная оценка плотности"),
         
         lty = 1:2, col = c("blue","red"), cex = .7)
  
}

residplot(OP3)

#Независимость остатков OP3

durbinWatsonTest(OP3)

#Линейность OP3

crPlots(OP3)

#Гомоскедастичность OP3

ncvTest(OP3)

spreadLevelPlot(OP3)

#Общая проверка выполнения требований OP3

gvmodelOP3 <- gvlma(OP3)

summary(gvmodelOP3)

#Мультиколлинеарность OP3

vif(OP3)

sqrt(vif(OP3)) > 2 # Проблема?

#######Сравнение моделей OP2 и OP3#######

anova(OP1, OP2, OP3)

AIC(OP1, OP2, OP3)

stepAIC(OP1, direction = "backward")

stepAIC(OP2, direction = "backward")

stepAIC(OP3, direction = "backward")

IZM1 <- list(durbinWatsonTest = durbinWatsonTest(OP1)$p, ncvTest = ncvTest(OP1)$p, vif = t(sqrt(vif(OP1)) > 2))

IZM2 <- list(durbinWatsonTest = durbinWatsonTest(OP2)$p, ncvTest = ncvTest(OP2)$p, vif = t(sqrt(vif(OP2)) > 2))

IZM3 <- list(durbinWatsonTest = durbinWatsonTest(OP3)$p, ncvTest = ncvTest(OP3)$p, vif = t(sqrt(vif(OP3)) > 2))

vv <- rbind.data.frame(OP1 = t(IZM1), OP2 = t(IZM2), OP3 = t(IZM3))

AIC <- AIC(OP1,OP2,OP3)$AIC

vv <- cbind.data.frame(vv, AIC)

Globalpvalue <- data.frame(pvalueOP1 = gvmodelOP1$GlobalTest$GlobalStat4$pvalue,
                           
                           pvalueOP2 = gvmodelOP2$GlobalTest$GlobalStat4$pvalue,
                           
                           pvalueOP3 = gvmodelOP3$GlobalTest$GlobalStat4$pvalue)

vv <- cbind.data.frame(vv, t(Globalpvalue))

#########################################


#savehistory()
save.image()