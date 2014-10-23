source("~/GitHub/SmartGeo/functions.R")

set.seed(123)
n_s = 100
n_t = 10
sp = data.frame( x=rnorm(n_s, sd=1) )
sp$y = ifelse(rbinom(n_s, size=1, p=.5)==0, sp$x, 100-sp$x)
ggplot(sp, aes(x=x, y=y) ) + geom_point()
t = data.frame(t=as.Date("2000-01-01") + 0:(n_t-1))

#Fit unidirectional variogram
d = merge(sp, t )
d$z = rnorm(n_s*n_t)
coordinates(d) = c("x", "y")
vst.s.only = variogram(z ~ 1, data=d, alpha=c(0,90))

temp = data.frame(sp)
max(dist(temp[temp$y>50,c("x","y")])); max(dist(temp[temp$y<50,c("x","y")]))

coordinates(sp) = c("x","y")
stdf = STFDF(sp=sp, time=t[,1], data=data.frame(d)[,"z",drop=F] )
vst.dir.old = gstat::variogramST( formula=z ~ 1, data=stdf, alpha=c(45,135), boundaries=0:20/2, cutoff=10)
vst.dir.new = variogramST( formula=z ~ 1, data=stdf, alpha=c(45,135), boundaries=0:20/2, cutoff=10)

vst.dir.manual1 = gstat::variogramST(z ~ 1, data=stdf[stdf@sp$y<50,], boundaries=0:20/2, cutoff=10)
vst.dir.manual2 = gstat::variogramST(z ~ 1, data=stdf[stdf@sp$y>50,], boundaries=0:20/2, cutoff=10)

p.s = ggplot( vst.s.only, aes(x=dist, y=gamma) ) + geom_line()
p.old = ggplot( vst.dir.old, aes(x=dist, y=gamma, color=timelag, group=timelag) ) + geom_line()
p.new1 = ggplot( vst.dir.new[!is.na(vst.dir.new$dir.hor) & vst.dir.new$dir.hor==45,]
    , aes(x=dist, y=gamma, color=timelag, group=timelag) ) + geom_line() + guides(color=F)
p.new2 = ggplot( vst.dir.new[!is.na(vst.dir.new$dir.hor) & vst.dir.new$dir.hor==135,]
    , aes(x=dist, y=gamma, color=timelag, group=timelag) ) + geom_line() + guides(color=F)
p.man1 = ggplot( vst.dir.manual1, aes(x=dist, y=gamma, color=timelag, group=timelag) ) +
  geom_line() + guides(color=F)
p.man2 = ggplot( vst.dir.manual2, aes(x=dist, y=gamma, color=timelag, group=timelag) ) +
  geom_line() + guides(color=F)

add = coord_cartesian(ylim=c(0,1.2), xlim=c(-.1,2.1))

grid.arrange( p.new1+add, p.new2+add, p.man1+add, p.man2+add, ncol=2)


################################################################
# Compare vst.dir.manual1/2 with vst.dir.new
################################################################

vst.dir.new = vst.dir.new[!is.na(vst.dir.new$np),]
vst.dir.new = vst.dir.new[vst.dir.new$np!=0,]
vst.dir.manual1 = vst.dir.manual1[!is.na(vst.dir.manual1$np),]
vst.dir.manual1 = vst.dir.manual1[vst.dir.manual1$np!=0,]
vst.dir.manual2 = vst.dir.manual2[!is.na(vst.dir.manual2$np),]
vst.dir.manual2 = vst.dir.manual2[vst.dir.manual2$np!=0,]

sum( vst.dir.manual1$np[vst.dir.manual1$timelag==0] )
n_s_low = sum( sp$y<50 )
n_s_low*(n_s_low-1)/2*n_t
sum(vst.dir.new$np[vst.dir.new$dir.hor==45 & vst.dir.new$timelag==0], na.rm=T)
vst.dir.new[vst.dir.new$dir.hor==45 & vst.dir.new$timelag==0,1:3] - 
  vst.dir.manual1[vst.dir.manual1$timelag==0,1:3]

sum( vst.dir.manual2$np[vst.dir.manual2$timelag==0] )
n_s_hi = sum( sp$y>50 )
n_s_hi*(n_s_hi-1)/2*n_t
sum(vst.dir.new$np[vst.dir.new$dir.hor==135 & vst.dir.new$timelag==0], na.rm=T)
vst.dir.new[vst.dir.new$dir.hor==135 & vst.dir.new$timelag==0,1:3] - 
  vst.dir.manual2[vst.dir.manual2$timelag==0,1:3]


sum( vst.dir.manual1$np[vst.dir.manual1$timelag==1] )
n_s_lo = sum( sp$y<50 )
#Now any station at current time and any station at lag time form a unique pair.
#But, only n_t-1 times
n_s_lo*n_s_lo*(n_t-1)
sum(vst.dir.new$np[vst.dir.new$dir.hor==45 & vst.dir.new$timelag==1
                   & vst.dir.new$dist<50], na.rm=T)
#gamma at dist=0 will be off, as manual uses a subset but new uses all pairs.
vst.dir.new[vst.dir.new$dir.hor==45 & vst.dir.new$timelag==1 & vst.dir.new$dist < 50,1:3] -
vst.dir.manual1[vst.dir.manual1$timelag==1 & vst.dir.manual1$dist < 50,1:3]

sum( vst.dir.manual2$np[vst.dir.manual2$timelag==1] )
n_s_hi = sum( sp$y>50 )
#Now any station at current time and any station at lag time form a unique pair.
#But, only n_t-1 times
n_s_hi*n_s_hi*(n_t-1)
#Slightly higher, as it has all pairs at dist=0
sum(vst.dir.new$np[vst.dir.new$dir.hor==135 & vst.dir.new$timelag==1
                   & vst.dir.new$dist<50], na.rm=T)
vst.dir.new[vst.dir.new$dir.hor==135 & vst.dir.new$timelag==1 & vst.dir.new$dist < 50,1:3] -
vst.dir.manual2[vst.dir.manual2$timelag==1 & vst.dir.manual2$dist < 50,1:3]
