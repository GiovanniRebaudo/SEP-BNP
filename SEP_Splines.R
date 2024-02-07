# Plot B-splines
# Codes accompanying "Separate Exchangeability as 
# Modeling Principle in Bayesian Nonparametrics"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(splines)
library(ggplot2)
library(latex2exp)
library(reshape2)
theme_set(theme_bw(base_size = 14))

x            <- seq(0, 1, by=0.001)
knots_seq    <- seq(0, 1, by=0.2)
oreder_spl   <- 4
spline_basis <- bs(x, degree=oreder_spl-1, knots=knots_seq)

# Simple plot in R
plot(spline_basis[,1]~x, ylim=c(0,max(spline_basis)), type='l', lwd=2, col=1, 
     xlab="Cubic B-spline basis", ylab="")
for (j in 2:ncol(spline_basis)){
  lines(spline_basis[,j]~x, lwd=2, col=j)
}

# Plot in R with ggplot2
plot_df <- melt(data.frame(x=x, spline_basis), id.vars="x", 
                variable.name="basis", value.name="y")
P = ggplot(plot_df) +
      aes(x=x, y=y, color=basis, group=basis) +
      geom_line() + ylab(latex2exp::TeX("$B_{i,4}(x)$")) +
      scale_color_discrete(guide=FALSE)

if(FALSE){
  ggsave(plot=P, file ="Image/Basis_Spline.pdf", 
         width=20, height=10, units = 'cm')
} else{
  P
}

# We have the following number of cubic splines
length(knots_seq) + 2

