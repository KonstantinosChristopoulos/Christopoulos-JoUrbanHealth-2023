#Read README file
#set up libraries#
library(rethinking)

#load data (not applicable)
D<-read.csv(".../gun.csv")

#transform data
#standardize exposure and covariates
D$Gstd<-standardize(D$gun)
D$PCstd<-standardize(D$property)
D$VCstd<-standardize(D$violent)

#create population offet
D$log_pop<-log(D$pop)

#create rates for the graphs
D$rate<-(D$deaths*100000)/D$pop

#create index variable for states
D$id<-coerce_index(D$state)

#create the mean gun onwership for each state
gbar<-sapply(1:50, function(j) mean(D$gun[D$id==j]))
D$Gbar<-gbar

#add to a list
d<-list(
	D=D$deaths,
	G=D$Gstd,
	Ge=D$Gsestd,
	VC=D$VCstd,
	PC=D$PCstd,
	state=D$state,
	P=D$log_pop,
	Gbar=D$Gbar,
	id=D$id
	)

#Mundlak models
#no other covariates (refered as "unadjusted" in the paper)
mmG<-ulam(
	alist(
		D~dgampois(lambda,phi),
		log(lambda)<- P + a[id] + b1*G + b5*Gbar[id],
		a[id]~dnorm(abar,sigma),
		abar~dnorm(-10,1),
		sigma~dexp(1),
		c(b1,b5)~dnorm(0,1),
		phi~dexp(1)
	),data=d, chains=4, cores=4 ,control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)
	
#crime rate adjsuted	
mmF<-ulam(
	alist(
		D~dgampois(lambda,phi),
		log(lambda)<- P + a[id] + b1*G + b2*VC + b3*PC + b4*VC*PC+ b5*Gbar[id],
		a[id]~dnorm(abar,sigma),
		abar~dnorm(-10,1),
		sigma~dexp(1),
		c(b1,b2,b3,b4,b5)~dnorm(0,1),
		phi~dexp(1)
	),data=d, chains=4, cores=4 , control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

#fixed-effects models	
#unadjusted
fencG<-ulam(
	alist(
		D~dgampois(lambda, phi),
		log(lambda)<- P+ a[id] +b1*G,
		b1~dnorm(0,1),
		a[id]~dnorm(-10,1), 
		phi~dexp(1) 
	),data=d, chains=4, cores=4 ,iter=4000, warmup=1000, log_lik=TRUE)	
	
#crime adjusted
fencF<-ulam(
	alist(
		D~dgampois(lambda, phi),
		log(lambda)<- P+ a[id] + b1*G + b2*VC + b3*PC +b4*PC*VC,
		a[id]~dnorm(-10,1), 
		c(b1,b2,b3,b4)~dnorm(0,1),
		phi~dexp(1) 
	),data=d, chains=4, cores=4 , iter=4000, warmup=1000, log_lik=TRUE)
	
post_F_fe<-extract.samples(fencF)
dens(post_F_fe$b1)
text(x,y,cex=0.8, "text",col="black")


#Figure 2a	
#unadjusted posterior of b1
post_G_fe<-extract.samples(fencG) 
post_G_mm<-extract.samples(mmG)

dens(post_G_mm$b1)
dens(post_G_fe$b1, col=rangi2, add=TRUE)
abline(v=0, lty=2)

#Figure 2c
#crime adjusted posterior of b1
post_F_fe<-extract.samples(fencF) 
post_F_mm<-extract.samples(mmF)

dens(post_F_mm$b1)
dens(post_F_fe$b1, col=rangi2, add=TRUE)
abline(v=0, lty=2)

#Figures 2b, 2d
#Posterior predictve simulation
#calculate mean posterior for the constant
mean(post_G_fe$a)
mean(post_F_fe$a)
mean(post_G_mm$a)
mean(post_F_mm$a)

#Figure 2b
#posterior predictive unadjusted
post_G_fe<-extract.samples(fencG) 
post_G_mm<-extract.samples(mmG)
lambdaMM<-function(G) exp(-10.42401 + post_G_mm$b1*G) 
lambdaFE<-function(G) exp(-10.50757 + post_G_fe$b1*G)
ns<-100
G.seq<-seq(from=-3, to=3 , length.out=ns)

#Mundlak unadjusted
muMM<-sapply(G.seq,lambdaMM)
lmuMM<-apply(muMM,2,mean)
lciMM<-apply(muMM,2,PI)
lmu1kMM<-lmuMM*100000
lci1kMM<-lciMM*100000
#FE unadjusted
muFE<-sapply(G.seq,lambdaFE)
lmuFE<-apply(muFE,2,mean)
lciFE<-apply(muFE,2,PI)
lmu1kFE<-lmuFE*100000
lci1kFE<-lciFE*100000

###plot
plot(NULL, xlab="Gun ownership (std)", ylab="Firearm deaths per 100k", 
	 xlim=c(-3,3) ,ylim=c(0,12))
#Mundlak
lines(G.seq, lmu1kMM ,lty=1 ,lwd=1.5)
shade(lci1kMM, G.seq ,xpd=TRUE)
#FE
lines(G.seq, lmu1kFE ,lty=2 ,lwd=1.5)
shade(lci1kFE, G.seq ,xpd=TRUE)

#Figure 2d
#posterior predictive adjusted
post_F_fe<-extract.samples(fencG) 
post_F_mm<-extract.samples(mmG)
lambdaMM<-function(G) exp( -10.54276 + post_F_mm$b1*G) 
lambdaFE<-function(G) exp(-10.49441 + post_F_fe$b1*G)
ns<-100
G.seq<-seq(from=-3, to=3 , length.out=ns)
#Mundlak
muMM<-sapply(G.seq,lambdaMM)
lmuMM<-apply(muMM,2,mean)
lciMM<-apply(muMM,2,PI)
lmu1kMM<-lmuMM*100000
lci1kMM<-lciMM*100000
#FE
muFE<-sapply(G.seq,lambdaFE)
lmuFE<-apply(muFE,2,mean)
lciFE<-apply(muFE,2,PI)
lmu1kFE<-lmuFE*100000
lci1kFE<-lciFE*100000

###plot
plot(NULL, xlab="Gun ownership (std)", ylab="Firearm deaths per 100k", 
	 xlim=c(-3,3) ,ylim=c(0,12))
#MM
lines(G.seq, lmu1kMM ,lty=1 ,lwd=1.5)
shade(lci1kMM, G.seq ,xpd=TRUE)
#FE
lines(G.seq, lmu1kFE ,lty=2 ,lwd=1.5)
shade(lci1kFE, G.seq ,xpd=TRUE)

#Figure 3 forest plot
#extract posteriors
forest<-list(
	MM_unadjusted=post_G_mm$b1,
	FE_unadjusted=post_G_fe$b1,
	MM_adjusted=post_F_mm$b1,
	FE_adjusted=post_F_fe$b1
	)

plot(precis(forest))

	
