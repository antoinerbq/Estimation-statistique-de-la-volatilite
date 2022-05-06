library(hawkes)
rm(list=ls())
set.seed(1)

#Parametres du hawkes:
lambda<-0.016
alpha<-0.024
beta<-0.11

#Simulation de Hawkes dependants:
l<-c(lambda,lambda)
a<-matrix(c(0,alpha,alpha,0),byrow=TRUE,nrow=2)
b<-c(beta,beta)
h<-3600*24 #Sur une journee.
Hawkes<-simulateHawkes(l,a,b,h)

#Hawkes contient les temps de sauts des deux processus de hawkes,
#on les extrait dans deux vecteurs pour reconstruire le processus.

t1<-Hawkes[[1]]
t2<-Hawkes[[2]]

#On cree donc une fonction qui renvoie la valeur au temps t d'un processus de comptage.
#En entree: le vecteur des temps de sauts, et le temps t.

val<-function(x,t){return(max((match((min(x[x>=t*(t<=x[length(x)])])),x)-1),length(x)*(t>x[length(x)])))}

#Commencons par verifier que le package ne renvoie pas quelquechose d'absurde,
#pour cela on teste si en prenant alpha et beta nuls on retrouve bien une Poisson(lambda*t) avec t=523.

z<-seq(0,1,0.01)
testP<-c()
for( i in 1:1000){
  hori<-3600
  x4<-simulateHawkes(lambda,0,0,hori)[[1]]
  testP<-c(testP,val(x4,523))
}

#On fait un diagramme quantile-quantile:

quantile(testP,probs = seq(0,1,0.1))
qqplot(testP,qpois(z,lambda*523),type="b")

x1<-seq(0,50,0.01)
lines(x1,x1,col="red")

#Celui-ci semble concluant, on décide de continuer d'utiliser ce package.

#Nous pouvons donc maintenant obtenir nos deux processus de Hawkes N1 et N2,
#et ainsi obtenir le processus de prix en faisant leur difference.
#Nous allons en creer deux: le premier prix sera observe toutes les 10 secondes sur 24h
#Le second sera observe toutes les 180 secondes sur 24h

ti<-seq(0,h,10)
gi<-seq(0,h,3*60)

N1_ti<-sapply(ti,val,x=t1)
N2_ti<-sapply(ti,val,x=t2)

N1_gi<-sapply(gi,val,x=t1)
N2_gi<-sapply(gi,val,x=t2)

prix_ti<-225+(N1_ti-N2_ti)
prix_gi<-225+(N1_gi-N2_gi)

#Nous avons debute les prix en 225 pour eviter d'avoir des prix negatifs
#car ceci ammene des difficultes pour le calcul des estimateurs.

#Observation de nos prix:


split.screen(1:2)
screen(1);plot(ti,prix_ti,type="s",xlab="t (secondes), une journee ecart 10 sec",ylab="X(t)")
screen(2);plot(gi,prix_gi,type="s",xlab="t (secondes),une journee ecart 3 min",ylab="X(t)")



#Observation de nos prix sur une heure afin d'observer de plus pres la "mean reversion":
zi<-seq(0,3600,10)
N1_zi=sapply(zi,val,x=t1)
N2_zi=sapply(zi,val,x=t2)
prix_zi=225+(N1_zi-N2_zi)


split.screen(1:2)
screen(1);plot(zi,prix_zi,type="s",xlab="t (secondes), la premiere heure",ylab="X(t)")
screen(2);plot(ti,prix_ti,type="s",xlab="t (secondes) sur une journee",ylab="X(t)")


#Maintenant que nous avons cree nos prix, nous pouvons en estimer la volatilite.
#Commencons par le plus simple: la volatilite realisee.
#On cree une fonction qui renvoie la valeur de la volatilite realisee selon la frequence d'observation des prix.
#En arguments: t qui represente l'horizon de temps sur lequel on observe un prix,
#tau qui est la frequence a laquelle on observe le prix en secondes.
rea_vol<-function(tau,t){
  wi<-seq(1,t,tau)
  N1_wi<-sapply(wi,val,x=t1)
  N2_wi<-sapply(wi,val,x=t2)
  prix_wi<-225+(N1_wi-N2_wi)
  return((1/t)*(sum(diff(prix_wi)^2)))
}

#Observations de la valeur de la volatilite realisee pour differents tau.
rea_vol(10,h)
rea_vol(120,h)
rea_vol(120*2,h)
rea_vol(120*3,h)

#On peut deja apercevoir une decroissance de la valeur suivant tau.
#Faisons donc un signature plot  (cf memoire partie 6 pour une definition),
#pour voir si nous arrivons a retrouver le premier paradoxe evoque.

tau<-seq(1,200,1)
sign<-sapply(tau,rea_vol,t=h)
plot(sign,type="l",xlab=expression(tau(s)),ylab="Signature plot")

#Rajoutons la courbe representant l'allure que devrait avoir le signature plot.
#La fonction a ete obtenue dans le papier:
#Modeling microstructure noise with mutually exciting points processes,
#De M.Hoffmann, E.Bacry, S.Delattre, J-F.Muzy .

C<-function(tau,a,b,mu){
  l<-((2*mu)/(1-(a/b)))
  k<-1/(1+(a/b))
  return(l*(k^2+(1-k^2)*((1-exp(-(a+b)*tau))/((a+b)*tau))))
}

curve(C(x,alpha,beta,lambda),add=TRUE,type="s",col="red",from=0,to=200)

#L'allure correspond bien, les resultats obtenus par simulation semblent donc coherents.

#Passons a l'estimateur de Fourier.
#Comme explique, pour pouvoir utiliser l'estimateur de Fourier nous devons
#passer l'intervalle de temps sur [0,2pi] .

FTime1<-((2*pi)/h)*ti
FTime2<-((2*pi)/h)*gi


#Premier estimateur de la volatilite integree par Fourier.
#Arguments: les log prix P,
#l'intervalle de temps reduit a [0,2pi] sur lequel on observe les prix,
#Le nombre de coefficients de Fourier a calculer pour reconstruire l'estimateur.
#Comme decrit dans le memoire, cet algorithme utilise le fait que:
#c_{-s}(dX)=Conj(c_{s}(dX)) rendant ainsi la somme symetrique.

int_vol<-function(P,t,N){
  rd<-diff(P)
  coeff0<-sum(rd)
  coeffp<-numeric(N)
  for(k in 1:N){
    coeffp[k]<-sum(exp(-1i*k*t[1:length(t)-1])*rd)
  }
  return(Re((1/(2*N+1))*(coeff0*Conj(coeff0)+2*sum(coeffp*Conj(coeffp)))))
  
}
#Commencons par essayer de determiner le choix optimal de N.
#Pour cela, nous faisons un graphique representant la valeur de l'estimateur par rapport a N:

si<-seq(1,500,5)
VolN<-sapply(si,int_vol,P=log(prix_ti),t=FTime1)
plot(VolN)

#On observe un plateau quand N devient grand.
#Nous choisissons donc une valeur de N assez grande pour avoir une valeur stable.
#Nous avons alors choisi de prendre N=500.

#On va regarder cette fois la valeur de l'estimateur de fourier comme une fonction de tau,
#qui est l'ecart de temps entre les prix.
four_vol<-function(tau,freq){
  wi<-seq(1,3600*24,tau)
  C1<-sapply(wi,val,x=t1)
  C2<-sapply(wi,val,x=t2)
  t<-((2*pi)/h)*wi
  P<-225+(C1-C2)
  r<-diff(log(P))
  c0<-sum(r)
  cp<-numeric(freq)
  for(k in 1:freq){
    cp[k]<-sum(exp(-1i*k*t[1:length(t)-1])*r)
  }
  return(Re((1/(2*freq+1))*(c0*Conj(c0)+2*sum(cp*Conj(cp)))))
}


r<-seq(1,200,1)

Vol1<-sapply(r,four_vol,freq=500)
plot(Vol1,type="l",ylim=c(0.03,0.045),xlab="expression(tau)(s)",ylab="Estimation de la volatilitÃ© par Fourier")

#On remarque que l'estimateur est plus stable en petite echelle que celui de la volatilite realisee.

#Valeurs de l'estimateur selon nos deux processus de prix.
int_vol(log(prix_ti),FTime1,500)
int_vol(log(prix_gi),FTime2,500)

#Comparaison par rapport aux valeurs de la volatilite realisee.
rea_vol(10,h)
rea_vol(300,h)

#On voit qu elles sont proches a haute frequence mais assez differentes a plus basse frequence.

#Implementons l'estimateur de Fourier utilisant le noyau de Fejer.
#Celui-ci a les memes arguments que l'estimateur precedent.

int_vol_fejer<-function(P,t,N){
  rd<-diff(P)
  coeffp<-numeric(2*N+1)
  coeffconj<-numeric(2*N+1)
  estim<-numeric(2*N+1)
  for(k in -N:N){
    coeffp[k+N+1]<-sum(exp(-1i*k*t[1:length(t)-1])*rd)
    coeffconj[k+N+1]<-sum(exp(1i*k*t[1:length(t)-1])*rd)
    estim[k+N+1]<-(1-abs(k)/N)*coeffp[k+N+1]*coeffconj[k+N+1]
  }
  return(Re(sum(estim)/(N+1)))
}

#Valeurs de l'estimateur suivant nos prix.
int_vol_fejer(log(prix_ti),FTime1,500)
int_vol_fejer(log(prix_gi),FTime2,500)

#Les valeurs sont tres similaires a celles trouvees avec l'autre estimateur.