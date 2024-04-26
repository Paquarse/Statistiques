# PARTIE 1

# L'objet de cette premiere partie est de programmer en R le calcul recursif de la moyenne et de la
# variance de chaque variable explicative, dans le but de centrer et reduire les observations dans l'etape
# suivante.

# Nettoyage de l'espace de travail

rm(list = ls())

# initialisation

p <- 49 # Nombre de variables explicatives

# Vecteur des mpyennes
moyx <- matrix(nrow = 1, ncol = p, 0)

# Vecteur des variances
varx <- matrix(nrow = 1, ncol = p, 0)

n <- 0

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne

nfile <- file("covtype_app.csv", open = "r")

while(length(oneLine <- readLines(nfile, n = 1, warn = FALSE)) > 0) {
  
  zz <- scan(text = oneLine, quiet = TRUE, sep = ";")
  
  # Selectionner les variables explicatives, soit les variables de 1 à p
  X <- zz[1:p]
  
  # Calcul itératif de la moyenne et de la variance pour les variables explicatives
  n <- n + 1
  
  # Variance
  varx <- varx + (1/n) * (((n - 1)/n) * (X - moyx) * (X - moyx) - varx)
  
  # Moyenne
  moyx <- ((n - 1) * moyx + X) / n
  
}

close(nfile)

# Calcul du vecteur des écart-types

sigmax <- sqrt(varx)
round(sqrt(sum(sigmax^2)), 3)

# On obient la norme euclidienne du vecteur des écarts-types de 2095.148
# Toutefois, pour s'assurer du bon fonctionnement de notre boucle, on a 
# importer la base de données, covtype_app et calculer la moyenne afin de
# comparer les résultats obtenus. 

cover_app = read.csv("covtype_app.csv", sep = ";")
cover_app_X = cover_app[,1:p]
Moyenne = as.numeric(apply(cover_app_X, MARGIN = 2, FUN = mean))
# Calcul de la distance entre les vraies moyennes et celles obtenue récursivement
dist = sum((Moyenne - moyx)^2)
dist

# On peut remarquer une distance très faible. Cela signifie que les moyennes obtenues 
# récursivement sont très proches de celles obtenues sur les vraies données d'apprentissage
# cover_type. On comprend donc la convergence de cet algorithmes vers les vraies valeurs
 
rm(cover_app, cover_app_X, Moyenne, dist)

# PARTIE 2

# L'objet de cette partie est de mettre en oeuvre l'algorithme de Newton stochastique pour 
# estimer recursivement les parametres du modele de regression logistique.

# Definition de la fonction pi(z)

pi <- function(x){
  return(1/(1+exp(-x)))
}

# Initialisation

p = 49

lambda = 1

Q <- lambda*diag(p+1)

theta <- rep(0,p+1)

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne

nfile = file("covtype_app.csv", open = "r")

while(length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz <- scan(text=oneLine, quiet=T, sep=";")
  
  # Variables explicatives
  
  X = zz[1:p]
  
  # Variable cible
  
  Y = zz[p+1]
  
  # centrage-réduction des variables explicatives 
  
  phi = c(1, as.numeric((X-moyx)/sigmax))

  # Mise à jour récursive de theta (Avec les nomenclature du cours)
  
  An = pi(as.numeric(t(theta)%*%phi))*(1-pi(as.numeric(t(theta)%*%phi)))
  
  Q = Q - An/as.numeric(1 + An*t(phi) %*% Q %*%phi) * (Q %*% phi %*% t(phi) %*% Q) 
  
  theta = theta + Q %*% phi*as.numeric(Y-pi(sum(t(theta)%*%phi)))
  
}

close(nfile)

# Norme du vecteur theta

round(sqrt(sum(theta^2)),3)

round(sum(sqrt(diag(Q))), 3)


# Calcul des intervalles de confiance à 95%

IC.binf = theta - qnorm(1 - 0.05/2, mean = 0) * sqrt(diag(Q)) # Vecteur des bornes inférieures

IC.bsup = theta + qnorm(1 - 0.05/2, mean = 0) * sqrt(diag(Q)) # Vecteur des bornes supérieures

# On retire la constante

IC.binf = IC.binf[-1]

IC.bsup = IC.bsup[-1]

# Identification des variables significatives

varsignif = (1:p)[0 < IC.binf | 0 > IC.bsup]

varsignif

sum(varsignif)

# Nous obtenons comme variables significatifs, les variables 1 à 10, 14, 15, 16, 19, 20, 31, 44 et 45
# Pour un total de 259. 

# PARTIE 3

# L'objet de cette partie est d'evaluer, sur l'echantillon test, les performances du modele
#  ajusté. Pour cela, nous utiliserons la table de confusion 

# initialisation

seuil = 0.5

N00 = 0

N01 = 0

N10 = 0

N11 = 0

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne

nfile = file("covtype_tst.csv", open = "r")

while (length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz = scan(text=oneLine, quiet=T, sep=";")
  
  # centrage-réduction des variables explicatives
  
  X = zz[1:p]
  
  Y = zz[p+1]
  
  phi =  c(1, as.numeric((X-moyx)/sigmax))
  
  # calcul de la prédiction
  
  ychap = t(theta) %*% phi
  
  if (pi(ychap) > seuil){
    if (Y == 0) {
      N01 = N01 + 1
    }
    if (Y == 1){
      N11 = N11 + 1
    }
  }
  
  if (pi(ychap) < seuil){
    if (Y == 0) {
      N00 = N00 + 1
    }
    if (Y == 1){
      N10 = N10 + 1
    }
  }
  
}

close(nfile)

N = N00 + N10 + N01 + N11

taux_erreurs = (N01+N10)/N

specificite = N11/(N10+N11)

sensibilite = N00/(N00+N01)


# Indices de performances

round(taux_erreurs*100,2)

round(specificite*100,2)

round(sensibilite*100,2)

# Matrice de confusion

confusion <- matrix(c(N00, N01, N10, N11), nrow = 2, byrow = TRUE,
                        dimnames = list(Cover_type_observe = c("0", "1"),
                                        Cover_type_predit = c("0", "1")))
t(confusion)


# Résumé des performances

resultat <- matrix(ncol=3,nrow=1)

colnames(resultat)<-c("Taux erreurs", "Specificités", "Sensibilité")

resultat[1,] = c(taux_erreurs,specificite,sensibilite)

resultat

# On a un taux d'erreur de 22%. Notre algorithme de newton stochastique donne 
# des résultats plus que correct. 

# Travail complémentaire

# Partie centrer que variables quantitatives

# Initialisation

p = 49

lambda = 1

Q<-lambda*diag(p+1)

theta <- rep(0,p+1)

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne
# D'apres l'enoncé, seule les 10 premieres sont numeriques, les autres sont binaires


nfile = file("covtype_app.csv", open = "r")

while(length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz <- scan(text=oneLine, quiet=T, sep=";")
  
  # Variables quantitatives : 1 à 11
  X_quant = zz[1:10]
  
  # Variable binaires
  X_bin = zz[11:p]
  
  # Variable cible
  Y = zz[p+1]
  
  #On centre et réduit les variables quantitatives avec la moyenne et écart type calculés itérativement
  
  X_quant_norm = as.numeric((X_quant-moyx[1:10])/sigmax[1:10])
  
  phi = c(X_quant_norm,X_bin)
  
  phi = c(1, phi)
  
  alpha <- pi(as.numeric(t(theta)%*%phi))*(1-pi(as.numeric(t(theta)%*%phi)))
  
  Q <- Q - alpha/as.numeric(1 + alpha*t(phi)%*% Q %*%phi) * (Q %*% phi %*% t(phi) %*% Q)
  
  theta <- theta + Q%*% phi*as.numeric(Y-pi(sum(t(theta)%*%phi)))
  
}



close(nfile)

# initialisation

seuil = 0.5

N00 = 0

N01 = 0

N10 = 0

N11 = 0

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne

nfile = file("covtype_tst.csv", open = "r")

while (length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz = scan(text=oneLine, quiet=T, sep=";")
  
  # centrage-réduction des variables explicatives
  
  X = zz[1:p]
  
  Y = zz[p+1]
  
  phi = phi = c(1, as.numeric((X-moyx)/sigmax))
  
  # calcul de la prédiction
  
  ychap = phi %*% theta
  
  if (pi(ychap) >= seuil){
    if (Y == 0) {
      N01 = N01 + 1
    }
    else{
      N11 = N11 + 1
    }
  }
  
  if (pi(ychap) <= seuil){
    if (Y == 0) {
      N00 = N00 + 1
    }
    else{
      N10 = N10 + 1
    }
  }
  
}

close(nfile)

N = N00 + N10 + N01 + N11

taux_erreurs = (N01+N10)/N

specificite = N11/(N10+N11)

sensibilite = N00/(N00+N01)

resultat2 <- matrix(ncol=3,nrow=1)

colnames(resultat2)<-c("Taux erreurs", "Specificités", "Sensibilité")

resultat2[1,] = c(taux_erreurs,specificite,sensibilite)

resultat
resultat2

# On remarque une légère augmentation du taux d'erreur quand on centre et réduit 
# les variables. La spécificité qui mésure la bonne prédiction des couvertures de 
# type sapin diminue. La sensibilité quant à elle augmente très très légèrement. 
# En conclusion, on peut remarquer le processus de centrage-réduction des variables 
# n'améliore pas les indicateurs. Ainsi, avec l'algorithme de Newton stochastique,
# il serait plus judicieux de ne pas centrer-reduire les variables. 

# Algorithme du gradient stochastique

#Initialisation

p = 49

lambda = 1

theta <- rep(0,p+1)

n = 0

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne

nfile = file("covtype_app.csv", open = "r")

while(length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz <- scan(text=oneLine, quiet=T, sep=";")
  
  X = zz[1:p]
  
  Y = zz[p+1]
  
  n = n + 1
  
  #On centre et réduit les données avec la moyenne et écart type calculés itérativement
  
  phi = c(1, as.numeric((X-moyx)/sigmax))
  
  gamma = n ^(-0.9)
  
  theta = theta - gamma * phi %*% (pi(t(theta) %*% phi) - Y)
  
}

close(nfile)

# initialisation

seuil = 0.5

N00 = 0

N01 = 0

N10 = 0

N11 = 0

# Connexion au fichier puis boucle de lecture séquentielle, ligne par ligne
nfile = file("covtype_tst.csv", open = "r")

while (length(oneLine <- readLines(nfile, n=1, warn=F)) > 0) {
  
  zz = scan(text=oneLine, quiet=T, sep=";")
  
  # centrage-réduction des variables explicatives
  
  X = zz[1:p]
  
  Y = zz[p+1]
  
  phi = phi = c(1, as.numeric((X-moyx)/sigmax))
  
  # calcul de la prédiction
  
  ychap = phi %*% theta
  
  if (pi(ychap) >= seuil){
    if (Y == 0) {
      N01 = N01 + 1
    }
    else{
      N11 = N11 + 1
    }
  }
  
  if (pi(ychap) <= seuil){
    if (Y == 0) {
      N00 = N00 + 1
    }
    else{
      N10 = N10 + 1
    }
  }
  
}
close(nfile)

N = N00 + N10 + N01 + N11

taux_erreurs = (N01+N10)/N

specificite = N11/(N10+N11)

sensibilite = N00/(N00+N01)

resultat3 <- matrix(ncol=3,nrow=1)

colnames(resultat3)<-c("Taux erreurs", "Specificités", "Sensibilité")

resultat3[1,] = c(taux_erreurs,specificite,sensibilite)

resultat3

# Avec l'algorrithme de gradien stochastique, nous obtenons à nouveau un taux d'erreur de 22% comme l'algorithme 
# de newton stochastique sans centré et réduire les données. La spécificité a augmenté est est maintenant à 75% 
# et la sensibilité est passée de 81% à 78%. 

# Dans ce cas précis, l'utilisation de cet algorithme n'améliore pas nettement nos résultats. D'autres valeurs de 
# gamma pourrait peut être améliorer ses résultats. 


# Pour aller plus loin, on a décidé de reprendre le même proccessus afin d'obtenir les variables significatifs

IC.binf = theta - qnorm(1 - 0.05/2, mean = 0) * sqrt(diag(Q)) # Vecteur des bornes inférieures

IC.bsup = theta + qnorm(1 - 0.05/2, mean = 0) * sqrt(diag(Q)) # Vecteur des bornes supérieures

# On retire la constante

IC.binf = IC.binf[-1]

IC.bsup = IC.bsup[-1]

# Identification des variables significatives

varsignif = (1:p)[0 < IC.binf | 0 > IC.bsup]

varsignif

sum(varsignif)

# Avec l'algorithme du gradient stochastique, on obtient comme variables explicatives, les
# variables 1  2  3  4  5  6  7  8 10 24. Cela représente certes, beaucoup moins de variables
# mais avec des performances presque aussi bonne que l'algorithme de newton stochastique. 
# le choix de l'algorithme de newton stochastique, pourrait donc se réveler, être plus parcimonieux
# Fin
