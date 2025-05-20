# data calibration curve
x <- rep(c(1.77, 3.54, 7.08, 14.16, 21.25, 28.33), each=5)
y <- c(5102, 4256, 5280, 4471, 4686, 13520, 10672,
       9470, 9801, 7600, 18564, 19577, 19637, 19913,
       21431, 53163, 53771, 46988, 47584, 62780, 88342,
       83381, 86648, 89513, 94384, 95578, 106670, 113805,
       102022, 123272)
data01 <- data.frame(cbind(x,y))

# weighted linear regression
lm1 <- lm(y~x, weights = 1/x^2)
plot(x,y)
abline(lm1, lty=2)
coeff1 <- as.data.frame(coef(lm1))
intercept <- coeff1[1,1]
slope <- coeff1[2,1]

#back-calculation 
x_pred<- (y - intercept)/slope
percent.bias <- (x_pred - x)/x*100
plot(x, percent.bias)
abline(h=0, lty=2)

# data subsetting 
res1 <- as.data.frame(cbind(x, percent.bias))
res1
res1.lloq <- subset(res1, x==1.77, percent.bias)
res1.3.54 <- subset(res1, x==3.54, percent.bias)
res1.7.08 <- subset(res1, x==7.08, percent.bias)
res1.14.16 <- subset(res1, x==14.16, percent.bias)
res1.21.25 <- subset(res1, x==21.25, percent.bias)
res1.28.33 <- subset(res1, x==28.33, percent.bias)

# calculation %-error
percent.error.lloq <- mean(res1.lloq$percent.bias)
percent.error.3.54 <- mean(res1.3.54$percent.bias)
percent.error.7.08 <- mean(res1.7.08$percent.bias)
percent.error.14.16 <- mean(res1.14.16$percent.bias)
percent.error.21.25 <- mean(res1.21.25$percent.bias)
percent.error.28.33 <- mean(res1.28.33$percent.bias)
mean.percent.error <- as.data.frame(rbind(percent.error.lloq,
                                          percent.error.3.54, percent.error.7.08,
                                          percent.error.14.16, percent.error.21.25,
                                          percent.error.28.33))
colnames(mean.percent.error)<-"Mean %-error"
mean.percent.error

# calculation rsd
conc1 <- as.data.frame(cbind(x, x_pred))
conc1
conc1.lloq <- subset(conc1, x==1.77, x_pred)
conc1.3.54 <- subset(conc1, x==3.54, x_pred)
conc1.7.08 <- subset(conc1, x==7.08, x_pred)
conc1.14.16 <- subset(conc1, x==14.16, x_pred)
conc1.21.25 <- subset(conc1, x==21.25, x_pred)
conc1.28.33 <- subset(conc1, x==28.33, x_pred)

rsd1.lloq <- sd(conc1.lloq$x_pred)/mean(conc1.lloq$x_pred)*100
rsd1.3.54 <- sd(conc1.3.54$x_pred)/mean(conc1.3.54$x_pred)*100
rsd1.7.08 <- sd(conc1.7.08$x_pred)/mean(conc1.7.08$x_pred)*100
rsd1.14.16 <- sd(conc1.14.16$x_pred)/mean(conc1.14.16$x_pred)*100
rsd1.21.25 <- sd(conc1.21.25$x_pred)/mean(conc1.21.25$x_pred)*100
rsd1.28.33 <- sd(conc1.28.33$x_pred)/mean(conc1.28.33$x_pred)*100

rsd1 <- as.data.frame(rbind(
  rsd1.lloq, 
  rsd1.3.54, 
  rsd1.7.08, 
  rsd1.14.16, 
  rsd1.21.25, 
  rsd1.28.33
))
colnames(rsd1) <- "%RSD"
rsd1

# data visualization

# calculation mean and sd
summary_data <- aggregate(y ~ x, data = data01, 
                          FUN = function(z) c(mean = mean(z), sd = sd(z),  
                                              n = length(z)))

# Memisahkan kolom rata-rata dan standar deviasi
summary_data$mean_y <- summary_data$y[, "mean"]
summary_data$sd_y <- summary_data$y[, "sd"]
summary_data$n_y <- summary_data$y[, "n"]
summary_data$se_y <- summary_data$sd_y / sqrt(summary_data$n_y)


# calculation CI from linear model
confidence_level <- 0.95
alpha <- 1 - confidence_level
t_critical <- qt(1 - alpha/2, df = summary_data$n_y - 1)
summary_data$ci_lower <- summary_data$mean_y - t_critical * summary_data$se_y
summary_data$ci_upper <- summary_data$mean_y + t_critical * summary_data$se_y

newx <- as.data.frame(x)

newy <- predict(lm1, interval = "confidence", level=.95, newdata = newx)
data02 <- data.frame(cbind(x, newy[,-1]))
y.lwr <- aggregate(lwr ~ x, data = data02, FUN = function(z) c(mean = mean(z), sd = sd(z),  n = length(z)))
y.upr <- aggregate(upr ~ x, data = data02, FUN = function(z) c(mean = mean(z), sd = sd(z),  n = length(z)))

# visualization with ggplot2
library(ggplot2)

p1 <- ggplot(summary_data, aes(x = x, y = mean_y)) +
  geom_point(size = 3) +  # Displaying the average data point.
  geom_errorbar(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y)) + # add error bar (sd)
  geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], linetype = "dashed", color = "blue") + # add regression line
  geom_line(data = y.lwr, aes(x = x, y = lwr[,1]), linetype = "dotted", color = "red")+
  geom_line(data = y.upr, aes(x = x, y = upr[,1]), linetype = "dotted", color = "red")+
  labs(
    title = "calibration curve of levofloxacin",
    x = "Levofloxacin in plasma (mcg/mL)",
    y = "Area Under Curve (Mean ± SD)"
  ) +
  theme_bw()
p1
