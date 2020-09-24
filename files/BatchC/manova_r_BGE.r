library(dplyr)
data <- read.csv('/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/manova_r_BGE_select.csv')

data$concentration<-factor(data$concentration, levels = c(1e-7, 1e-5, 1e-3, 1e-1), labels = c("CBT", "BT", "NT", "HT"))
str(data)


dataMelt <- melt(data, id=c('concentration'), measured=c('BGE', 'r'))



names(dataMelt) <- c('Concentration', 'Response', 'Frequency')

ocdBoxplot <- ggplot(dataMelt, aes(Concentration, Frequency, color = Response))
ocdBoxplot + geom_boxplot() + labs(x='Concentration', y='', color='Response') + scale_y_log10()


outcome <- cbind(data$BGE, data$r)
ocdModel <- manova(outcome ~ concentration, data=data)
summary(ocdModel, intercept=TRUE)
