#Lecture des données

rm(list=ls())

# Question 1
data = read.table(file="Boston_housing.data")

n = length(data[,1])
p = ncol(data)-1

colnames(data)[1:p] = paste("X",1:p,sep="")
colnames(data)[(p+1)] = "Y"

# Estimation avec la fonction lm
reslm <- lm(Y ~.,data)
summary(reslm)
coef.lm = summary(reslm)$coeff[,1]
s2.lm = (summary(reslm)$sigma)^2

# Méthode récursive

delta = 1e6
Q = delta * diag(p+1)
theta = rep(0, (p+1))
s2 = 0

Q <- delta * diag(p + 1)


for (j in 1:n) {
  
  y <- data[j, p + 1]
  
  phi <- c(1, as.numeric(data[j, 1:p]))
  
  # Mise à jour de theta et s2
  
  U <- Q %*% phi
  
  f <- as.numeric(t(phi) %*% U)
  
  R <- y - as.numeric(t(theta) %*% phi)
  
  s2 <- s2 + (1/j)*((1/(1 + f))*(R%*%R)  - s2)
    
  Q <- Q -  (as.numeric(1/(1 + t(phi) %*% U)) * (U %*% t(U)))
  
  theta <- theta + Q %*% phi * (y - as.numeric(t(phi) %*% theta))
  
}


# Calcul des erreurs d'estimation
erreur.coef = sqrt(sum((coef.lm - theta)^2))
erreur.coef
erreur.sigma2 = abs(s2.lm - s2*n/(n-p-1))
erreur.sigma2




# Comparaison des coefficients estimés
Result = round(cbind(coef.lm, theta) , 5)
colnames(Result) = c("MC ordinaires","MC récursifs")
round( Result, 4)

xtable::xtable(round(Result, 4))



##  Question  2

estimate_errors <- function(delta_values, data) {
  n <- nrow(data)
  p <- ncol(data) - 1
  
  errors <- matrix(NA, nrow = length(delta_values), ncol = 3,
                   dimnames = list(NULL, c("Delta", "Erreur.coefficients", "Erreur.sigma2")))
  
  for (k in seq_along(delta_values)) {
    delta <- delta_values[k]
    Q <- delta * diag(p + 1)
    theta <- rep(0, p + 1)
    s2 <- 0
    
    for (j in 1:n) {
      y <- data[j, p + 1]
      phi <- c(1, as.numeric(data[j, 1:p]))
      
      # Mise à jour de theta et s2
      U <- Q %*% phi
      f <- as.numeric(t(phi) %*% U)
      R <- y - as.numeric(t(theta) %*% phi)
      s2 <- s2 + (1 / j) * ((1 / (1 + f)) * (R %*% R) - s2)
      Q <- Q - (as.numeric(1 / (1 + t(phi) %*% U)) * (U %*% t(U)))
      theta <- theta + Q %*% phi * (y - as.numeric(t(phi) %*% theta))
    }
    
    # Calcul des erreurs d'estimation
    erreur.coef <- sqrt(sum((coef.lm - theta)^2))
    erreur.sigma2 <- s2 * n / (n - p - 1)
    
    errors[k, ] <- c(delta, erreur.coef, erreur.sigma2)
  }
  
  return(errors)
}

# Exemple d'utilisation avec différentes valeurs de delta
delta_values <- c(1e-5, 0.2, 10, 100, 1000, 10000, 1e6, 1e8, 1e10, 1e11)  # Différentes valeurs de delta à tester
errors <- estimate_errors(delta_values, data)
errors

xtable::xtable(errors)

library(ggplot2)

# Convertir les données en un data frame
errors_df <- as.data.frame(errors)

# Tracer les courbes
ggplot(errors_df, aes(x = Delta)) +
  geom_line(aes(y = Erreur.coefficients, color = "Erreur.coefficients")) +
  geom_line(aes(y = Erreur.sigma2, color = "Erreur.sigma2")) +
  scale_x_log10() +  
  labs(x = "Delta", y = "Erreur", color = "Type d'erreur") +
  theme_minimal()


# Partie 2 

# Nettoyage de l'espace de travail

close(nfile)
rm(list=ls())

# Ouverture du Fichier de données

nfile  = file("data_MCR.txt", open = "r")

# Lecture de la première ligne de données

oneLine = readLines(nfile, n = 1,  warn = FALSE)

x = scan(text=oneLine, quiet=T, sep=";")

p = length(x) - 1 # Nombre de régresseurs


# Initialisation

delta =  1e5 # a remplir

Q      = diag(p+1) * delta

Beta = rep(0, (p+1))

# Initialisation de l'estimateur récursif  de la variance sigma2
s2 = 0

# Initialisation des matrices utiles
Z = rep(0, p+1)

S = matrix(0, ncol = p+1, nrow=p+1)

Xbar = rep(0,p)

SX = matrix(0, ncol=p, nrow=p)

nval = 1

while (nval < 20000) {
  
  # construction des données
  y   = x[p+1]
  
  x = x[1:p]
  
  phi = c(1, x)
  
  # Mise à jour de Z et de S
  
  S = S + phi %*% t(phi)
  
  Z = Z + phi*y
  
  # Mise à jour de  Xbar, SX
  
  SX = SX + ((nval-1)/nval)*tcrossprod(x- Xbar)
  
  Xbar = Xbar + (x-Xbar)/nval
  
  # Mise à jour de beta et s2
  
  U <- Q %*% phi
  
  fn = as.numeric( t(phi) %*% Q %*% phi )
  
  R <- y - as.numeric(t(Beta) %*% phi)
  
  s2 = s2 + (1/nval)*((1/(1 + fn))*(R%*%R)  - s2)
  
    Q = Q - (as.numeric(1/(1 + t(phi) %*% U)) * (U %*% t(U)))
    
    Beta = Beta + Q %*% phi * (y - as.numeric(t(phi) %*% Beta))
    
    # Acquisition d'une nouvelle donnée
  
      oneLine = readLines(nfile, n = 1,  warn = FALSE)
  
  x = scan(text=oneLine, quiet=T, sep=";")
  
  nval = nval + 1

}

# Calcul de l'estimateur direct

B = solve(S) %*% Z

# Comparaison des deux estimations
sqrt(sum((Beta-B)^2))



# valeur a mettre dans le CR

alpha = 0.05

# Calcul des intervalles de confiance à 95%

IC.binf = Beta - qnorm((1 - alpha/2), mean = 0, sd = 1)*sqrt(diag(Q))

IC.bsup = Beta + qnorm((1 - alpha/2), mean = 0, sd = 1)*sqrt(diag(Q))


# On retire la constante
IC.binf = IC.binf[-1]
IC.bsup = IC.bsup[-1]

# Identification des variables significatives
varsignif = (1:p)[IC.binf > 0 | IC.bsup  < 0]
varsignif
# Valeurs a mettre dans le CR
sum(varsignif)
length(varsignif)

close(nfile)

# 4 - Calcul avec les données restantes

rm(list=ls())

nfile  = file("data_MCR.txt", open = "r")

# Lecture de la première ligne de données

oneLine = readLines(nfile, n = 1,  warn = FALSE)

x = scan(text=oneLine, quiet=T, sep=";")

p = length(x) - 1 # Nombre de régresseurs

# Initialisation des variables pour EQM et S2Y

EQM = 0

S2Y = 0

y_mean = 0

delta = 0.5  
Q = diag(p+1) * delta
beta = rep(0, (p+1))

nval = 20001
n = 1

while(TRUE) {
  oneLine = readLines(nfile, n = 1, warn = FALSE)
  if(length(oneLine) == 0) { break }  # Quitte la boucle si aucune ligne n'est lue
  x = scan(text=oneLine, quiet=TRUE, sep=";")
  if(length(x) != p + 1) { next }  # Passe cet
  
  
  # Extraction de la variable à expliquer (y) et des régresseurs (x)
  y = x[p+1]
  x = x[1:p]
  phi = c(1, x)
  
  # Calcul de y_mean
  y_mean = y_mean + (1/n) * (y - y_mean)
  
  # Calcul de la variance à expliquer (S2Y)
  
  S2Y = S2Y + ((y - y_mean)^2)
  
  
  fn <-t(phi) %*% Q %*% phi
  Q_phi <- Q %*% phi
  Q_phi_t_phi_Q <- Q_phi %*% t(phi) %*% Q
  r<-1 + fn
  r<-as.numeric(r)
  term1 <- Q_phi_t_phi_Q * (1/r)
  Q<- Q - term1
  pred <-as.numeric(t(phi) %*% beta)
  
  # Calculer l'erreur de prédiction
  rn <- y - pred
  
  # Mettre à jour beta
  beta <- beta + Q %*% phi * (rn)
  
  # Calcul de l'erreur quadratique moyenne (EQM)
  
  EQM = EQM + (rn^2)
  
  n = n + 1
  
}

EQM=EQM/80000
S2Y=S2Y/80000

close(nfile)

EQM
S2Y

# Calcul du pourcentage de variance expliquée

pve = (1 - (EQM / S2Y)) * 100
pve















