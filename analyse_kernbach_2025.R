################################################################################
######################## Oekologie der Lebensraeume 2025 #######################
################################################################################

# loescht alles bisherige aus der Global Environment
rm(list = ls(all = TRUE))
install.packages("vegan")
install.packages("ggplot2")

# Liste an Paketen definieren, die geladen werden sollen
packages <- c("vegan", "ggplot2")

# die Funktion library fuer die gesamte Liste Packages ausfuehren
lapply(packages, library, character.only = TRUE)

# setzt Pfad zum Verzeichnis in dem die Daten liegen
workFolder <- normalizePath(getwd(), winslash="\\", mustWork=TRUE)
setwd(workFolder)

# Einlesen der Daten
kernbach.arten <-read.csv2("Rohdaten Kernbach(Tabelle1).csv",header = TRUE) 

# Tabelle anzeigen
head(kernbach.arten)
names(kernbach.arten)
fix(kernbach.arten)

# Entfernen erste Spalte (Taxonomie)
kernbach.arten <- kernbach.arten[,-(1)]

# Artname als "row names"
row.names(kernbach.arten) <- kernbach.arten$Art

# Löschen der Spalte "Artname"
kernbach.arten <- kernbach.arten[,-(1)]

kernbach.arten

# Anlegen als data frame
kernbach.arten.t <- as.data.frame(t(kernbach.arten)) 
d <- c("K", "F")
d <- rep.int(d,9)
kernbach.arten.t$kolk <- d

str(kernbach.arten.t)

# Sortierung der Daten
S_kernbach <- specnumber(kernbach.arten.t[,c(-39)])### Artenzahl
H_kernbach <- diversity(kernbach.arten.t[,c(-39)])### Shannon-Diversität
J_kernbach <- H_kernbach/log(specnumber(kernbach.arten.t[,c(-39)]))# Evenness
Abu_kernbach <- rowSums(kernbach.arten.t[,c(-39)])### Abundanz=Anzahl der Individuen pro Probestelle

# Anhängen der Indices an den data frame
kernbach.arten.t$S <- S_kernbach
kernbach.arten.t$H <- H_kernbach
kernbach.arten.t$J <- J_kernbach
kernbach.arten.t$Abu <- Abu_kernbach

kernbach.arten.t

# Erstellen der Spalte mit Enfernungen zur Quelle - hier einfach Reihenfolge
kernbach.arten.t$dist <- c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)

kernbach.arten.t


# Erstellung von Subsets für die fließende und die aufgestauten Bereiche Bereiche
kernbach.arten.t_F <- subset(kernbach.arten.t, kolk == "F")
kernbach.arten.t_K <- subset(kernbach.arten.t, kolk == "K")

#Alpha-Diversität

# Leeres Plot-Fenster vorbereiten
plot(kernbach.arten.t$dist, kernbach.arten.t$S,
     type = "n",
     main = "Alphadiversität entlang des Kernbaches",
     xlab = "Distanz zur Quelle",
     ylab = "Anzahl der Arten")

# Punkte für Fließbereich (F) einzeichnen
points(kernbach.arten.t_F$dist, kernbach.arten.t_F$S,
       col = "blue", pch = 16)

# Punkte für Kolkbereich (K) einzeichnen
points(kernbach.arten.t_K$dist, kernbach.arten.t_K$S,
       col = "red", pch = 17)

# Regressionslinie für alle Daten
model1 <- lm(S ~ dist, data = kernbach.arten.t)
abline(model1, col = "black", lwd = 2)

# Legende hinzufügen
legend("topright",
       legend = c("Fließ", "Kolk"),
       col = c("blue", "red"),
       pch = c(16, 17))

# Lade ggplot2
library(ggplot2)

# Erstelle das Regressionsmodell
model1 <- lm(S ~ dist, data = kernbach.arten.t)
class(model1)
anova(model1)

# Residuen extrahieren
residuen <- model1$residuals

# Modellzusammenfassung
summary(model1)

# Histogramm der Residuen mit Achsen und Farbe
hist(model1$residuals,
     main = "Verteilung der Residuen",
     xlab = "Residuen",
     col = "lightblue",
     border = "white",
     breaks = 10)

# Beta-Diversität

# Erstellen einer Matrix der Paarweisen Unterschiede der Artenzusammensetzung
kernbach.arten.t.beta <- kernbach.arten.t

#Entfernen nicht benötigter Spalten
kernbach.arten.t.beta <- kernbach.arten.t.beta[,-c(44:49)]

kernbach.arten.t.beta <- as.data.frame(kernbach.arten.t.beta)
kernbach.arten.t.beta <- na.omit(kernbach.arten.t.beta)
str(kernbach.arten.t.beta)

# nicht numerische Daten ausschließen
daten_nur_numeric <- kernbach.arten.t.beta[, sapply(kernbach.arten.t.beta, is.numeric)]

vegdist(daten_nur_numeric, binary = TRUE, method = "jaccard")

beta_div_kernbach <- vegdist(daten_nur_numeric, binary = TRUE, method = "jaccard")

# Jaccard Distanzen hier nur presence/absence

# - Erstellen einer Matrix der paarweisen Distanzen

# 1. Nur die Artendaten extrahieren (Presence/Absence)
arten_daten <- kernbach.arten.t[, sapply(kernbach.arten.t, is.numeric)]
arten_daten <- arten_daten[, 1:38]  

# 2. Jaccard-Distanzmatrix berechnen
library(vegan)
beta_div_kernbach <- vegdist(arten_daten, method = "jaccard", binary = TRUE)

# 3. Geographische Distanzmatrix (Entfernung zur Quelle)
# Prüfen, ob Spalte 'dist' existiert und numerisch ist
if (!"dist" %in% names(kernbach.arten.t)) {
  stop("Spalte 'dist' fehlt im Datensatz.")
}
if (!is.numeric(kernbach.arten.t$dist)) {
  stop("'dist' muss numerisch sein.")
}
geo_dist <- dist(kernbach.arten.t$dist)

# 4. Mantel-Test
mantel_result <- mantel(beta_div_kernbach, geo_dist, method = "pearson", permutations = 1000)

# 5. Ergebnis anzeigen
print(mantel_result)

# 6. Visualisierung: Jaccard-Distanz vs. geographische Distanz
# Umwandlung der Distanzobjekte in Vektoren für Plot
geo_dist_vec <- as.vector(geo_dist)
beta_div_vec <- as.vector(beta_div_kernbach)

plot(geo_dist_vec, beta_div_vec,
     main = "Zusammenhang zwischen geographischer Distanz und Beta-Diversität",
     xlab = "Geographische Distanz (zur Quelle)",
     ylab = "Jaccard-Distanz (Beta-Diversität)",
     pch = 1,
     col = "black")

# Lade ggplot2
library(ggplot2)

# Bereite die Daten vor
kernbach.arten.t.dist <- kernbach.arten.t

# Berechne die Distanzmatrix (z. B. geografische Distanz oder Flussdistanz)
dist_matrix <- as.matrix(dist(kernbach.arten.t.dist$dist))

# Berechne die Beta-Diversität (z. B. Jaccard-Distanz)
beta_matrix <- as.matrix(dist(kernbach.arten.t.dist$S))  # Falls S die Artenzusammensetzung ist

# Erstelle einen DataFrame mit allen Paaren
df <- data.frame(
  dist = as.vector(dist_matrix),
  beta = as.vector(beta_matrix)
)

# Lade ggplot2
library(ggplot2)

# Erstelle das Regressionsmodell
model1 <- lm(beta ~ dist, data = df)

# Residuen extrahieren
residuen <- model1$residuals

# Histogramm der Residuen mit Base R
hist(residuen,
     main = "Histogramm der Residuen",
     xlab = "Residuen",
     col = "lightblue",
     border = "white")

# Mantel's Test
mantel(beta_div_kernbach,geo_dist,  method = "pearson",
       permutations = 1000)

# Kombination aus Alpha- und Beta-Diversität

# Bereite die Daten vor
kernbach.arten.t.dist <- kernbach.arten.t

# Berechne Distanzmatrix (z. B. geografische Distanz)
dist_matrix <- as.matrix(dist(kernbach.arten.t.dist$dist))

# Berechne Beta-Diversität (z. B. Jaccard-Distanz)
beta_matrix <- as.matrix(dist(kernbach.arten.t.dist$S))  # S = Artenzusammensetzung

# Erstelle DataFrame mit allen Paaren
df <- data.frame(
  dist = as.vector(dist_matrix),
  beta = as.vector(beta_matrix)
)

# Füge Alpha-Diversität hinzu (z. B. Mittelwert der beiden Probestellen)
alpha_values <- kernbach.arten.t.dist$S  # Alpha-Diversität pro Standort
alpha_matrix <- outer(alpha_values, alpha_values, FUN = function(x, y) (x + y) / 2)
df$alpha <- as.vector(alpha_matrix)

# Modell 1: Nur Beta-Diversität
model_beta <- lm(dist ~ beta, data = df)

# Modell 2: Alpha + Beta-Diversität
model_combined <- lm(dist ~ alpha + beta, data = df)

# Modelle definieren
model_beta <- lm(dist ~ beta, data = df)
model_alpha <- lm(dist ~ alpha, data = df)
model_combined <- lm(dist ~ alpha + beta, data = df)

# Zusammenfassungen anzeigen
summary(model_beta)
summary(model_alpha)
summary(model_combined)

# Modellvergleich 1: Beta vs. Kombiniert
anova(model_beta, model_combined)

# Modellvergleich 2: Alpha vs. Kombiniert
anova(model_alpha, model_combined)