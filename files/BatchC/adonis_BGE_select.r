dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))
### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)

library(dplyr)
data <- read.csv('/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/adonis_BGE_select.csv')

X <- data.frame(data$response, as.factor(data$rX), data$concentration)
names(X) <- c('Response', 'rX', 'Concentration')
Y <- data.frame(data$BGE)
names(Y) <- c('BGE')
row_sub = apply(Y, 1, function(row) all(row !=0 ))
Yf <- Y[row_sub,]
Xf <- X[row_sub,]
mod <- metaMDS(Yf)
plot(mod)
with(Xf, ordihull(mod, group=rM, show="1"))
with(Xf, ordihull(mod, group=rM, show="2", col=3))
with(Xf, ordispider(mod, group=Response, lty=3, col="red"))

adonis(Yf ~ rX, data=Xf, perm=999)
