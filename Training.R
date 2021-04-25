list.dirs("D:/RData/DIPLOM")
library("dplyr")
library("caret")
library("AER")
library("ggplot2")
library("sandwich")
library("ivpack")

h <- read.csv("D:/RData/DIPLOM/dataFull.csv", header=TRUE, sep=";")
glimpse(h)

OP1 <- lm(VVP ~ SG4Z + X4BR, data = h)
OP2 <- lm(VVP ~ X4BRZ + SNZP, data = h)

in_train <- createDataPartition(y = h$VVP, p=0.75, list=FALSE) # для обучения берём 75% данных

# Для дальнейшей оценки модели методом МНК
h_train <- h[in_train,] # сюда берём только наблюдения для обучающей выборки
h_test <- h[-in_train,] # оценка качества прогнозов, обучающие данные исключаем, остальные оставляем

nrow(h)
nrow(h_train)
nrow(h_test)

model_1 <- lm(data=h_train, VVP ~ SG4Z + X4BR)
model_2 <- lm(data=h_train, VVP ~ X4BRZ + SNZP) # высокая точность

y <- h_test$VVP

# Прогнозируем
y_hat_1 <- predict(model_1, h_test)
y_hat_1 # Очень низкая точность прогноза
y_hat_2 <- predict(model_2, h_test)
y_hat_2 # Точность очень высокая. Для 2017 года реальное ВВП = 92101, а в прогнозе = 92878.077

# Сумма квадратов ошибок прогнозов
sum((y-y_hat_1)^2)
sum((y-y_hat_2)^2)


nextYear <- data.frame(X4BRZ=2000, SNZP=55000)
nextPredict <- predict(model_2, newdata=nextYear)
nextPredict
