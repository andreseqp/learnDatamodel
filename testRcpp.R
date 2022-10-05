library(Rcpp)
library(here)
library("BayesianTools")

sourceCpp("test.cpp")

timesTwo(3)

?VSEM

focal<-new(model_param)
focal$set_gamma(0.5,0.5)
focal$logist()
new_param<-new(model_param)
new_param$set_gamma(1,10)
new_param$logist()

focal$copy(new_param$get_ptr())

focal$logist()


focal$set_gamma(1,1)
focal$logist()
focal$set_gamma(0.5,0.6)
focal$logist()

hello()
bla()
bla2(42, 0.42)

w <- new(World)
w$greet()
w$set("hohoho")
w$greet()

muUnif<-new(Uniform,0,2)

muUnif$draw(5)

sourceCpp("ActCrit_R.cpp")
