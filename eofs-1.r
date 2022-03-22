#
# Empirical Orthgonal Functions (EOFs) mit SVD berechnen
#   - Vergleich zu Diagonalisierung mit Eigenvektoren
#   - Laufzeiten ~m^2*n fÃ¼r n>=m
#   - Daten-Kompression/Fehler
#   - Implementierung: DGESSD in Lapack
#

##
# Notation nach: 
# B. Halldor and Venegas, S.A., (1997).
# A manual for EOF and SVD analyses of climate data.
# McGill University, CCGCR Report No. 97-1, MontrÃ©al, QuÃ©bec, 52pp. 
##

#
# synthetisiere skalare Zeitreihen: 
#    nx * ny rÃ¤umliches Gitter, 
#    m Zeitpunkte
#

m=600    # LÃ¤nge der Zeitreihen (Statistik: "StichprobengrÃ¶ÃŸe")
nx=ny=20 # Gitter-Dimensionierung
F.3d=array(dim=c(m,nx,ny)) # Anlegen im Speicher
dim(F.3d)

# Daten-Synthese erfolgt durch sin/cos Terme und Rauschen
x=(1:nx)/nx*pi  # x-Werte (Gitter-Koordinaten)
y=(1:ny)/ny*pi  # y-Werte (Gitter-Koordinaten)

t=(1:m)/m*2*pi  # Zeit-Werte

for (i in 1:m) # im Kern der Schleife: ganze y-Vektoren 
 for (j in 1:nx)  F.3d[i,j,]=(cos(t[i])*sin(x[j])*cos(y)+
          +1.5*sin(t[i]*2+y)*cos(x[j]-y*2) + 
          - sin(t[i]*6)*cos(x[j]*8)
          + rnorm(ny)/25) # normalverteilter Fehler

# Speichere 2D-Gitter Variable in einem Vektor, je Zeit
#  Zeile=Gitter-Vektor, Spalte = Zeitreihe in einem Gitterpunkt
F=matrix(F.3d,m,nx*ny) # hintere 2 Dimensionen (Gitter) 
#                        in 2. Dim. von F packen
# zeitliche Mittelwerte in jedem Gitterpunkt abziehen 
#   -> nx*ny Zeitreihen der "Anomalien"="zentrierte Daten"
for(i in 1:(nx*ny)) F[,i]=F[,i]-mean(F[,i])

# Ansicht erstes Grid, farbcodiert rot=min, weiss=max
image(matrix(F[1,],nx,ny))

# fÃ¼r Moviemaker: mit Abspeichern als .png
for (i in 1:(m/10)*10) {
    #png(paste(i,"png",sep="."));
    image(matrix(F[i,],nx,ny),main=i);
    #dev.off()
  }

# SVD Berechnen, Rechenzeit anzeigen
system.time(svd(F)->S)
str(S)  # Struktur: zwei Matrizen u,v, ein Vektor d
# Plot formatieren
par(mfrow=c(2,1),mar=c(2,1,1,1),oma=c(2,2,0,0))
# Darstellung erste 20 SingulÃ¤rwerte
plot(S$d[1:20],type="b",main="erste 20 Singulaer-Werte")
# Darstellung als Anteil an Gesamt-Varianz 
Var=S$d^2/sum(S$d^2)
plot(Var[1:20]*100,type="b",main="% Anteil an Gesamtvarianz")
#welche Singulaerwerte erscheinen nicht vernachlaessigbar?

M=1:5
# Anzeige der ersten M EOFs
par(mfrow=c(3,2)) # in einer Abb.
for (j in M){
#  png(paste("EOF.",j,".png"))
 image(matrix(S$v[,j],nx,ny),main=paste("EOF",j,":",round(100*Var[j],1),"% Var"))
 contour(matrix(S$v[,j],nx,ny),add=TRUE,labcex=1)
 #dev.off()
} 

# SpÃ¤ter im Vergleich mit Eigenvektoren der Cov.Matrix

# absoluter Fehler aus der Rueckrechnung 
error <- F - S$u %*% diag(S$d) %*% t(S$v)
range(error)

# Zeitreihen (PCs)
M <- 1:5
par(mfrow=c(length(M),1),mar=c(2,1,1,2))
for (i in M) plot(S$u[,i],type="l",main=paste("PC ",i))


# Fehler bei Kompression auf nur 2 EOFs
M=1:2
error.2=F - S$u[,M] %*% diag(S$d[M]) %*% t(S$v[,M])
range(error.2)
range(F)

M=1:4
error.4=F - S$u[,M] %*% diag(S$d[M]) %*% t(S$v[,M])

# Verbindung zum Eigenwert-Problem
# Kovarianz-Matrix Z; bis auf Faktor n; siehe Vorlesung
t(F) %*% F-> Z
dim(Z)
#image(Z)
# Eigenwerte und Vektoren berechnen
system.time(eigen(Z)->E)
str(E)
# Gleichheit d_i^2 == lambda_i
range(E$values-S$d^2)
dim(E$vectors)
range(E$values)

# Multiplot-Abbildung formatieren
M=1:4
par(mfrow=c(length(M),2),mar=c(1,1,2,1),oma=c(1,1,1,0))

# Vergleich der EOFs
#  Eigenvektoren nur bis auf Vorzeichen bestimmt
for (i in M) {
  # Projektion erstes Gitter auf die Eigenvektoren e_M
  image(matrix(E$vectors[,i],nx,ny),main=paste(i,":",round(E$values[i],6)),axes=F)
  contour(matrix(E$vectors[,i],nx,ny),add=TRUE)
  image(matrix(S$v[ ,i],nx,ny),main=paste(i,":",round(S$d[i]^2,6)),axes=F)
  contour(matrix(S$v[ ,i],nx,ny),add=TRUE)
}

par(mfrow=c(4,2),mar=c(2,1,2,1),oma=c(1,1,2,0))
for (i in 1:4){
  plot(F%*%E$vectors[,i],main=paste("PC",i,"(Eigen)"),pch=".")
  plot(S$u[,i]*S$d[i],main=paste("PC",i,"(SVD)"),pch=".")}  

## Rechenzeit-Verhalten svd() (Lapack) in R
for(m in c(100,200,400,800,1600)){
  n=1600
  print(system.time(svd(matrix(rnorm(m*p),m,p))))}
# ~m^2  n>=m

for(n in c(400,800,1600,3200,6400)){
  m=400
  print(system.time(svd(matrix(rnorm(m*n),m,n))))}
#~n ,   n>=m


# fÃ¼r Movie der komprimierten Rekonstruktion, nur jedes 10. Bild
M=1:4
S$u[,M] %*% diag(S$d[M]) %*% t(S$v[,M])->F.red
for (i in 1:(m/10)*10) {
  #png(paste(i,"-compr-4EOFs.png"));
  image(matrix(F.red[i,],nx,ny));
  #dev.off()
}

#absoluter Fehler
for (i in 1:(m/10)*10) {#png(paste(i,"-error-4EOFs.png"));
  image(matrix(F[i,],nx,ny)-matrix(F.red[i,],nx,ny));
  #dev.off()
}

### WeiterfÃ¼hrendes 

# Experimente mit dem Ergebnis
#
# ist z=cov(D*V^T) diagonal?

image((diag(S$d) %*% S$v)%*%t(diag(S$d) %*% S$v))
z=((diag(S$d) %*% S$v)%*%t(diag(S$d) %*% S$v))
diag(z)[1:10] # nur wenige groÃŸe EintrÃ¤ge auf Diagonale

# ohne S: sind alle Diagonal-Werte diag(Cov(v^T)) ~ 1
all(abs(diag((S$v)%*%t(S$v))-1.0)<1e-14)

# Diagonalisierung ? t(C) cov(F) C
(t(E$vectors)%*%t(F)%*%F%*%(E$vectors))[1:4,1:4]

# siehe auch Video
# maximale Fehler, alle m Samples
par(mar=c(3,4,2,2),mfrow=c(1,1))
plot(apply(abs(F-F.red),1,function(p) max(p)),ylab="max(error)",xlab="sample" )

# Karte maximaler Fehler
ij.max=apply(abs(F-F.red),1,function(p) which.max(p))
#
xx=(ij.max %% nx )/(nx-1);yy=floor(ij.max/ny)/(ny-1) 
plot(xx,yy,col=densCols(xx,yy),pch=19,xlab="x",ylab="y",cex=2)

image(matrix(apply(F,1,sd),nx,ny))
contour(matrix(apply(F,1,sd),nx,ny),add=TRUE,nlevels = 5,labcex = 0.9)
# hellblau:   geringster Fehler
# dunkelblau:    grÃ¶ÃŸter Fehler
points(xx,yy,col=densCols(xx,yy),pch=19,xlab="x",ylab="y",cex=2)



#
# (C) Stephan Frickenhaus, 2016
#