install.packages("prob")
library(prob)
library(ggplot2)

t = seq(-10,10,.1)
vals = sapply(t, function(t)cfunif(t, min=-1,max=1) )
qplot(t, Re(vals) )
qplot(t, Im(vals) )

vals = sapply(t, function(t)cfnorm(t) )
qplot(t, Re(vals) )
qplot(t, Im(vals) )

qplot( t, cos(t)*exp(-t^2/1000) )

x = seq(-5,5,.1)
k = 1
vals = sapply(x, function(x){
  r = integrate( function(t){Re(cos(t)*exp(-1i*t*x-t^2/k))}
    ,lower=-Inf, upper=Inf )$value
  i = integrate( function(t){Im(cos(t)*exp(-1i*t*x-t^2/k))}
    ,lower=-Inf, upper=Inf )$value
  return(c(r,i))
} )/(2*pi)
qplot( x, vals[1,] )
qplot( x, vals[2,] )
qplot( x, vals[1,] ) + geom_line( aes(y=dnorm(x,sd=1.8)) )



#From stackoverflow:
characteristic_function_to_density <- function(
  phi, # characteristic function; should be vectorized
  n,   # Number of points, ideally a power of 2
  a, b # Evaluate the density on [a,b]
) {
  i <- 0:(n-1)            # Indices
  dx <- (b-a)/n           # Step size, for the density
  x <- a + i * dx         # Grid, for the density
  dt <- 2*pi / ( n * dx ) # Step size, frequency space
  c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
  d <-  n/2 * dt          # (center the interval on zero)
  t <- c + i * dt         # Grid, frequency space
  phi_t <- phi(t)
  X <- exp( -(0+1i) * i * dt * a ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}

#My own:
#' @param density Density function, should be a vectorized function of x
#' @param t Grid of t values at which to evaluate the characteristic function
#' @param ... Additional parameters to pass to f
density_to_char <- function(f, t, ...){
  out = sapply(t, function(tVal){
    return(c(
      integrate( function(x){Re(exp(-1i*tVal*x)*f(x,...))}, -Inf, Inf)$value
     ,integrate( function(x){Im(exp(-1i*tVal*x)*f(x,...))}, -Inf, Inf)$value
    ))
  } )
  out = out[1,] + out[2,]*1i
  return(out)
}

d <- characteristic_function_to_density(
  function(t,mu=1,sigma=.5) 
    exp( (0+1i)*t*mu - sigma^2/2*t^2 ),
  2^8,
  -3, 3
)
plot(d$x, d$density, las=1)
curve(dnorm(x,1,.5), add=TRUE)

d <- characteristic_function_to_density(
  function(t)sapply(t, function(t){cfunif(t,min=0,max=1)}),
  2^10,
  -1, 2
)
plot(d$x, d$density, pch=16, cex=.5)

d <- characteristic_function_to_density(
  function(t)
    cos(t)*exp(-t^2/100),
  2^10,
  -5, 5
)
plot(d$x, d$density, pch=16, cex=.5)

t = seq(-40,40,.4)
phi.unif = density_to_char( dunif, t=t, min=-1 )
qplot(t, Re(phi.unif), geom="line", color="Real") +
  geom_line(aes(y=Im(phi.unif), color="Imaginary") )

phi.abs = density_to_char( function(x){
    ifelse(x > -1 & x < 1, abs(x), 0)
  }
  ,t=t )
qplot(t, Re(phi.abs), geom="line", color="Real") +
  geom_line(aes(y=Im(phi.abs), color="Imaginary") )

2*integrate(exp, 0, 1)$value
2*(exp(1)-1)
phi.exp = density_to_char( function(x){
    ifelse(x > -1 & x < 1, exp(abs(x))/(2*(exp(1)-1)), 0)
  }
  ,t=t )
qplot(t, Re(phi.exp), geom="line", color="Real") +
  geom_line(aes(y=Im(phi.exp), color="Imaginary") )

2*integrate(sqrt, 0, 1)$value
phi.sqrt = density_to_char( function(x){
    ifelse(x > -1 & x < 1, sqrt(abs(x))/(4/3), 0)
  }
  ,t=t )
qplot(t, Re(phi.sqrt), geom="line", color="Real") +
  geom_line(aes(y=Im(phi.sqrt), color="Imaginary") )

integrate(function(x){abs(x)^9}, -1, 1)
phi.pow = density_to_char( function(x){
    ifelse(abs(x) < 1, abs(x)^9/.2, 0)
  }
  ,t=t )
qplot(t, Re(phi.pow), geom="line", color="Real") +
  geom_line(aes(y=Im(phi.pow), color="Imaginary") )


qplot(t, Re(phi.sqrt), geom="line", color="Sqrt(x)") +
  geom_line(aes(y=Re(phi.exp), color="Exp(x)") ) +
  geom_line(aes(y=Re(phi.abs), color="Abs(x)") ) +
  geom_line(aes(y=Re(phi.pow), color="Abs(x)^9") )

qplot(t, 1-Re(phi.sqrt), geom="line", color="Sqrt(Abs(x)) on [-1,1]") +
  geom_line(aes(y=1-Re(phi.exp), color="Exp(Abs(x)) on [-1,1]") ) +
  geom_line(aes(y=1-Re(phi.abs), color="Abs(x) on [-1,1]") ) +
  geom_line(aes(y=1-Re(phi.pow), color="Abs(x)^9 on [-1,1]") ) +
  labs(color="Original Density", title="Variogram Model")

t = seq(-700,700,.1)
const = integrate(function(x){abs(x)^100}, -1, 1)$value
phi.pow2 = density_to_char( function(x){
    ifelse(abs(x) < 1, abs(x)^100/const, 0)
  }
  ,t=t )
ggsave("~/Personal Files/Pretty Plots/characteristic_function_for_density_x_to_100_on_minus1_to_1.png",
  qplot(t, 1-Re(phi.pow2), geom="line", color="Variogram") +
    guides(color=F)
)

#Or, avoid numerical issues altogether:
# 