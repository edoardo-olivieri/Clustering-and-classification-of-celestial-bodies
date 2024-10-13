library(tidyverse)
library(skimr)
library(GGally)
library(mclust)
library(gridExtra)
library(caret)
library(flexmix)
##########ATTENZIONE############
#Alcuni blocchi di codice impiegano molto tempo per girare a causa della grande 
#quantità di dati.


#Ripuliamo i dati dalle variabili non rilevanti e lo salviamo in una tibble 
star<- as_tibble(read_csv("C:/.../Verna Olivieri Portanova 3.csv"))
star<- star %>% select(u,g,r,i,z,redshift,class)

#Usiamo il comando skim() della library skimr che è simile al comando summary(),
#fornisce alcune statistiche di base, tra cui il numero di valori mancanti per 
#ogni singola variabile e i vari quantili
skim(star)
star$class <- as.factor(star$class)
#Notiamo che non ci sono valori mancanti e che le variabili u,g,z e redshift 
#presentano un valore minimo insolito, molto lontano dal primo quantile

#Diamo un'occhiata ai boxplot delle variabili
as_tibble(star %>% select(-class)) %>%
  gather(variabili, valori) %>%
  ggplot(aes(x=reorder(variabili, valori, FUN=median), y=valori, fill=variabili)) +
  geom_boxplot(show.legend=F, fill="cornflowerblue", color="darkblue", outlier.shape="x", outlier.size=3) +
  facet_wrap(~variabili, scales="free")+
  labs(title="Boxplots - Variabili") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())
#Notiamo che è presente un outlier con un valore anomalo e decidiamo di eliminarlo
star<-star %>% filter(g > -3000)

#Controlliamo nuovamento i boxplot delle variabili
as_tibble(star %>% select(-class)) %>%
  gather(variabili, valori) %>%
  ggplot(aes(x=reorder(variabili, valori, FUN=median), y=valori, fill=variabili)) +
  geom_boxplot(show.legend=F, fill="cornflowerblue", color="darkblue", outlier.shape="x", outlier.size=3) +
  facet_wrap(~variabili, scales="free")+
  labs(title="Boxplots - Variabili") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())

#Controlliamo la distribuzione delle variabili attraverso i corrispettivi 
#istogrammi con sovraimposta funzione di densità stimata
star %>% 
  select(-class) %>%
  gather(variabili, valori) %>% 
  ggplot(aes(x=valori, y=..density..)) +
  geom_histogram(fill="cornflowerblue", color="white", alpha=0.3) +
  geom_density(color="cornflowerblue")+
  facet_wrap(.~variabili, scales = "free")+
  labs(title="Distribuzione - Variabili") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())

#Per completezza, visto che siamo in possesso delle etichette vere, riguardiamo 
#le precedenti rappresentazioni grafiche suddividendo per gruppi

######IMPIEGA MOLTO TEMPO PER PRODURRE UN RISULTATO######
ggpairs(star[,-c(7)], aes(color = star$class, alpha=0.5)) +
  scale_fill_manual(values = c("cornflowerblue","lightblue","darkblue")) +
  scale_color_manual(values = c("cornflowerblue","lightblue","darkblue"))

#-------------------------------------------------------------------------------

#Attravero l'analisi delle componenti principali della matrice di varianze e 
#covarianza individuiamo le variabili più rilevanti e creiamo una nuova tibble 
#con solo quest'ultime per svolgere le successive analisi
pca = princomp(star%>%select(-class),cor=T)
sum((pca$sdev[1:2])^2)/6 #Le prime due componenti spiegano circa il 90% della variabilità
nome_var<-names(star)[apply(pca$loadings[,1:2], 2, function(x) which.max(x^2))]
star_new<-star%>%select(all_of(nome_var), class)

#Siccome conosciamo le etichette diamo un'occhiata allo scatterplot di r-redshift 
#diviso per gruppi
ggplot()+
  geom_point(data=star, aes(x=redshift, y=r, color=class),alpha=0.5, size = 0.6, shape = 16)+
  scale_color_manual(values = c("cornflowerblue","lightblue","darkblue"), 
                     name = "Legenda", labels = c("GALAXY","QSO",  "STAR"))+
  guides(color = guide_legend(override.aes= list(size=2))) +
  labs(title="Scatterplot -  r-redshift")

#Procediamo con il model based clustering, individuando il modello e il numero 
#di componenti ottimale in base all'ICL

######IMPIEGA MOLTO TEMPO PER PRODURRE UN RISULTATO######
out.mclustICL1 <- mclustICL(star_new%>%select(-class))
summary(out.mclustICL1)   #VVV4   VVI8   VVE8
plot(out.mclustICL1, legendArgs = list(x = "bottomright"))
abline(v=4, h=max(as.numeric(summary(out.mclustICL1))), lty=2)

#Il modello migliore risulta essere il VVV con 4 gruppi e procediamo 
#nell'implementarlo
out.mclust1<-Mclust(star_new%>%select(-class), G=4, modelNames="VVV")
summary(out.mclust1)
star_new %>%
  ggplot(aes(x=redshift, y=r, color=as.factor(out.mclust1$classification)))+
  geom_point(alpha=1, show.legend=F, size = 0.6)+
  scale_color_manual(values = c("purple","lightblue","purple4","darkblue"))+
  labs(title="Scatterplot - 4 cluster trovati")

#Controlliamo la validità del modello calcolandone l'entropia e confrontadola 
#con il valore soglia
(entropia<-out.mclust1$bic-out.mclust1$icl)   #11069.85
n<-nrow(star_new)
n*log(4)

#Sapendo che i gruppi sono 3 ci chiediamo come mai il numero ottimale sia 4,
#decidiamo quindi di osservare meglio gli istogrammi con sovraimposta funzione 
#di densità stimata delle due variabili
d1<-star_new%>%
  mutate(cluster=star$class)%>%
  ggplot(aes(x=r, y=..density..))+
  geom_histogram(fill="cornflowerblue", color="white", alpha=0.3)+
  geom_density(color="cornflowerblue")+
  facet_wrap(.~cluster, scales = "free")+
  labs(title="Distribuzione r") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())
d2<-star_new%>%
  mutate(cluster=star$class)%>%
  ggplot(aes(x=redshift, y=..density..))+
  geom_histogram(fill="cornflowerblue", color="white", alpha=0.3)+
  geom_density(color="cornflowerblue")+
  facet_wrap(.~cluster, scales = "free")+
  labs(title="Distribuzione redshift") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())

grid.arrange(d1,d2,nrow=2)
#Notiamo che nel gruppo GALAXY vi è un andamento bimodale per entrambe le 
#variabili

#Proviamo a rifare il model based clustering imponendo 3 gruppi
out.mclustICL2 <- mclustICL(star_new%>%select(-class), G=3)
summary(out.mclustICL2)   #VVV3   VVI3   VVE3
out.mclust2<-Mclust(star_new%>%select(-class), G=3, modelNames="VVV")
summary(out.mclust2)

#Calcoliamo il CER 
lab.true<-star$class
lab.cluster<-out.mclust2$classification
classError(lab.cluster, class=lab.true)

#Rappresentiamo le osservazioni misclassified
miss_class<-classError(lab.cluster, class=lab.true)$misclassified
star%>%ggplot()+
  geom_point(aes(x=redshift, y=r, color=as.factor(out.mclust2$classification)),alpha=0.5, size = 0.6)+
  scale_color_manual(values = c("lightblue","darkblue","cornflowerblue","red"), 
                     name = "Legenda", labels = c("QSO","STAR","GALAXY","ERROR"))+
  labs(title="Scatterplot - 3 cluster trovati con errori sovraimposti")+
  guides(color = guide_legend(override.aes= list(size=2))) +
  geom_point(data=star[miss_class,], aes(x=redshift, y=r, color="red"),alpha=0.08, size = 0.6)

#Calcoliamo l'ARI
adjustedRandIndex (lab.cluster , lab.true)   #0.7897513

#Calcoliamo la confusion matrix

lab.true<-as.factor(lab.true)
lab.cluster <- as.factor(lab.cluster)
levels(lab.cluster)<-c("QSO", "STAR", "GALAXY")
confusionMatrix(lab.cluster, lab.true)

#-------------------------------------------------------------------------------

#Procediamo effettuando un Mixtures of Experts Model 
#Individuiamo il numero di componenti ottimale sulla base dell'ICL

set.seed(123)

######IMPIEGA MOLTO TEMPO PER PRODURRE UN RISULTATO######
fit1<-stepFlexmix(r~redshift, data=star_new, nrep=5, k=1:5, verbose=T, drop=F, unique=F, concomitant=FLXPmultinom(~redshift))
plot(ICL(fit1))
points(which.min(ICL(fit1)), min(ICL(fit1)), col="red", pch=19)

#Individuiamo il modello con 3 componenti migliore sulla base dell'ICL
iclval<-Inf
itermax<-5
set.seed(123)

######IMPIEGA MOLTO TEMPO PER PRODURRE UN RISULTATO######
for(i in 1:itermax){
  final<-flexmix(r~redshift, data=star_new, k=3, concomitant=FLXPmultinom(~redshift))
  if(iclval>ICL(final)){
    iclval<-ICL(final)
    bestfinal<-final
  }
}
(bestfinal)  #24234    49986    25779
                                                                                 
#Analizziamo l'experts network del modello
asymp.inf<-refit(bestfinal)
summary(asymp.inf)
plot(asymp.inf)

#Analizziamo la gating network del modello
summary(bestfinal)
parameters(bestfinal ,which="concomitant")

#Analizziamo il rootogram delle probabilità a posteriori
plot(bestfinal)

#Calcoliamo CER del modello
labs<-bestfinal@cluster
classError(labs, class=star_new$class)    #0.1237212

#Rappresentiamo graficamente i cluster ottenuti
ggplot(data=star_new, mapping = aes(y=r, x=redshift, color=factor(labs)))+
  geom_point(alpha=0.2, size=0.6)+ 
  geom_smooth(method="lm", linewidth=1, se=F)+
  scale_color_manual(values = c("darkblue","cornflowerblue","lightblue"),
                     name="Legenda", labels = c("STAR","GALAXY","QSO"))

#-------------------------------------------------------------------------------

#Attraverso un diagramma a torta controlliamo se il dataset è bilanciato
tab<-as.data.frame(table(star$class)) #Creiamo una tabella di freq. asso. delle classi
colnames(tab)<-c("Class", "Freq")
rt<-rev(tab$Freq)#Salviamo le freq. asso. in ordine crescente
pos<-cumsum(rt)-rt/2 #Calcoliamo e salviamo la posizione da dare alle etichette
lbs<-paste(round(rt/sum(rt)*100,2), "%", " (", rt, ")", sep="")#Calcoliamo e 
#salviamo le percentuali insieme alle freq. asso. nelle etichette

tab%>%
  ggplot(aes(x=factor(1), y=Freq, fill=Class))+
  geom_col()+
  coord_polar(theta="y", direction=-1)+ #Trasforma il geom_col in un diagramma a torta
  theme_void()+
  scale_fill_manual(values = c("cornflowerblue", "lightblue","darkblue"))+
  geom_label(x=1.4, y=pos, aes(label=lbs), fill="lightyellow", size=4)#Aggiungiamo le etichette
  

#Siccome non risulta bilanciato appllico al dataset downSample() della library caret
star_bilanciato<-downSample(star_new, as.factor(star_new$class))%>%select(-Class)

#Procediamo con la MDA
set.seed(123)
lista<-list()
lista[["er"]]<-c(1)

######IMPIEGA MOLTO TEMPO PER PRODURRE UN RISULTATO######
for (i in 1:5){
  #Usiamo createDataPartition della library caret per dividere il dataset in due
  partizione<-createDataPartition(y=star_bilanciato$class, p=0.8, list=F) 
  training<-star_bilanciato[partizione,]
  testing<-star_bilanciato[-partizione,]
  mod<-MclustDA(training%>%select(-class), training$class)
  er<-cvMclustDA(mod, nfold=10)$ce #V-fold CV con V=10
  if (lista$er>er){#Se er è migliore del precedente i risultati vengono salvati in una lista
    lista$er<-er
    lista[["training"]]<-training
    lista[["testing"]]<-testing
    lista[["mod"]]<-mod
    lista[["predizione"]]<-predict(mod,testing%>%select(-class))
    lista[["MER"]]<-classError(lista$predizione$classification, class=testing$class)
  }
}
(lista$er)   #0.0504318
mod<-lista$mod
(MER<-lista$MER)   #0.0524789
predizione<-lista$predizione
confusionMatrix(predizione$classification, as.factor(testing$class))

#Rappresentiamo graficamente i classification boundaries
prec <- 150
x1 <- seq(-0.1, 7.1, length = prec)
x2 <- seq(9.81, 29.7, length = prec)
s <- expand.grid(x2, x1)
P <- predict.MclustDA(mod, s)$z
s$group <- factor(max.col(P))
testing <- testing %>% mutate(class = as.factor(class))
pastel<-0.7
ggplot() +
  geom_point(data = s, aes(x = Var2, y = Var1, color = group), size = 3) +
  geom_point(data = testing, aes(x = redshift, y = r, color = class), size = 0.7) +
  scale_color_manual(values = c(rgb(1, pastel, pastel), rgb(pastel, 1, pastel), rgb(pastel, pastel, 1), 2, 3, 4)) +
  labs(title = "Classification boundaries", y = "r", x = "redshift")+
  theme(legend.position = "none")
