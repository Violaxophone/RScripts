################################################################################
######################## Oekologie der Lebensraeume 2025 #######################
################################################################################

################################################################################
#### General setttings #########################################################
################################################################################
#loescht alles bisherige aus der Global Environment
rm(list = ls(all = TRUE))

user_profile <- Sys.getenv("USERPROFILE")
workFolder <- file.path(user_profile, "Dokumente", "Studium", "Biologie", "6. Semester", "Ökologie der Lebensräume", "Exkursion", "Kernbach")

#setzt Pfad zum Verzeichnis in dem die Daten liegen
file.exists(workFolder)

setwd(normalizePath(workFolder))
install.packages("vegan")

#Liste an Paketen definieren, die geladen werden sollen
packages <- c("vegan", "ggplot2")

###Alternativ
##library(vegan)
##library(ggplot2)

#die Funktion library fuer die gesamte Liste Packages ausfuehren
lapply(packages, library, character.only = TRUE)

################################################################################
# I.    Kernbach
################################################################################

# Fragen:
# 1. Gibt es einen Zusammenhang zwischen Artenreichtum und Entfernung der 
# Probestellen von der Quelle?
# 2. Unterscheiden sich stark fließenden Bereiche und die aufgestauten Bereiche 
# in Individuenreichtum und Artenreichtum?
# 3. Ist die Artenzusammensetzung unterschiedlicher, je weiter zwei Standorte
# voneinander entfernt sind (isolation by distance)?

# Auswertung:
# (1) LM mit Artenreichtum als abhängige Variable, Entfernung zur Quelle als 
# kontinuierliche Variable - wenn L + S getrennt analysiert, oder GLM mit 
# Artenreichtum als abhängige Variable, Entfernung von der Quelle 
# als kontinuierliche unabhängige Variable, fließend/aufgestaut als Faktor. 
# (3) Manteltest: Matrix Entfernung gegen Jaccard-Index oder ähnliches.


################################################################################
# 1. Artenreichtum - Alphadiversität - kann mit der Funktion specnumber aus dem
# Paket vegan erfasst werden. Dafür wird eine Artentabelle benötigt.

kernbach.arten <-read.csv2("Rohdaten Kernbach(Tabelle1).csv",header = TRUE) 
# Um die Daten genauer anzusehen,eignen sich z.B. die folgenden Funktionen - 
# head, names, fix.
head(kernbach.arten)
names(kernbach.arten)
fix(kernbach.arten)
# Beim fix Befehl muss das Fenster manuell geschlossen werden.
# Für die weiteren Analysen benötigen wir nicht mehr die Namen der Ordnungen; Spalte
# 1 kann entfernt werden. Mit eckigen Klammern kann einzelne Zeilen, Spalten, 
# oder Zellen zugegriffen werden [Zeile,Spalte]. D.h. wir löschen einfach die 
# erste Spalte; Doppelpunkt heisst bis; hier also von Spalte 1 bis Spalte
# 2.
kernbach.arten <- kernbach.arten[,-(1)]

### Im nächsten Schritt weißen wir Spalte 1 mit den Arten als rownames aus

row.names(kernbach.arten) <- kernbach.arten$Art

#### und löschen dann Spalte 1 mit den Artnamen

kernbach.arten <- kernbach.arten[,-(1)]

kernbach.arten

# Im Folgenden werden wir das Paket Vegan nutzen, um Diversitätsmaße wie Arten-
# reichtum zu berechnen - mit den Funktionen specnumber, diversity, betadiver. 
# Dazu muss die Tabelle transponiert werden - wofür es in R die Funktion "t" 
# gibt. 
# Da R beim Transponieren aus dem data.frame eine Matrix macht, 
# geben wir explizit vor, dass das ganze als Datenframe angelegt wird. Das Er-
# gebnis lässt sich bspw. durch den is(...) Befehl überprüfen.
kernbach.arten.t <- as.data.frame(t(kernbach.arten)) 
d <- c("K", "F")
d <- rep.int(d,9)
kernbach.arten.t$kolk <-d

str(kernbach.arten.t)

#####write.csv2(kernbach.arten.t, "kernbach_2025_t.csv")


# Berechnung der Artenzahl:
# Anwendung specnumber Funktion (aus vegan) auf Artentabelle. Erste Spalte
# mit Namen der Aufnahmeflächen muss ausgelassen werden. Mit Hilfe der eckigen
# Klammern können wir auf einen Teil des Datensatzes zugreifen - in diesem Fall
# sparen wir die erste Spalte aus. Zeichen vor dem Komma bezieht sich auf 
# Zeilen. Zeichen nach dem Komma - auf Spalten. Wir schreiben das Ergebnis 
# direkt in eine neue Spalte "a.div" in unseren Datenframe. 
S_kernbach <- specnumber(kernbach.arten.t[,c(-39)])### Artenzahl
H_kernbach<-diversity(kernbach.arten.t[,c(-39)])### Shannon-Diversität
J_kernbach<-H_kernbach/log(specnumber(kernbach.arten.t[,c(-39)]))# Eveness
Abu_kernbach<-rowSums(kernbach.arten.t[,c(-39)])### Abundanz=Anzahl der Individuen pro Probestelle

# Anhängen der Indices an den data frame
kernbach.arten.t$S <-S_kernbach
kernbach.arten.t$H <-H_kernbach
kernbach.arten.t$J <-J_kernbach
kernbach.arten.t$Abu <-Abu_kernbach

kernbach.arten.t

# Erstellen der Spalte mit Enfernungen zur Quelle - hier einfach Reihenfolge
kernbach.arten.t$dist <- c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)

kernbach.arten.t


# Erstellung von Subsets für die fließende und die aufgestauten Bereiche Bereiche
kernbach.arten.t_F <- subset(kernbach.arten.t, kolk == "F")
kernbach.arten.t_K  <- subset(kernbach.arten.t, kolk == "K")

# unterschiedliche Farben Habitate
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

# plotten der Distanzen zur Quelle gegen alpha Diversität
plot(kernbach.arten.t$S~kernbach.arten.t$dist,
     main = "Alphadiversität entlang des Kernbaches",
     xlab= "Distanz zur Quelle",
     ylab = "Anzahl der Arten")
model1<-lm(kernbach.arten.t$S~kernbach.arten.t$dist)
summary(model1)
abline(model1)
hist(model1$residuals)

# Lade ggplot2, falls noch nicht geschehen
library(ggplot2)

# Erstelle das Regressionsmodell
model1 <- lm(S ~ dist, data = kernbach.arten.t)
class(model1)
anova(model1)

# Residuen extrahieren
residuen <- model1$residuals

# verschönerte Darstellun
# Scatterplot mit optimierter Darstellung
plot(kernbach.arten.t$dist, kernbach.arten.t$S,
     main = "Alphadiversität entlang des Kernbaches",
     xlab = "Distanz zur Quelle (m)",
     ylab = "Anzahl der Arten",
     pch = 19, col = "steelblue", cex = 1.4)

# Regressionslinie hinzufügen
model1 <- lm(kernbach.arten.t$S ~ kernbach.arten.t$dist)
abline(model1, col = "firebrick", lwd = 2)

# Modellzusammenfassung
summary(model1)




### Call:
###  lm(formula = kernbach.arten.t$S ~ kernbach.arten.t$dist)
###
### Residuals:
### Min      1Q  Median      3Q     Max 
### -3.1818 -1.6879 -0.4621  1.3364  6.2818 
###
### Coefficients:
### Estimate Std. Error t value Pr(>|t|)   
### (Intercept)             3.6455     1.0259   3.553  0.00227 **
### kernbach.arten.t$dist   0.5121     0.1922   2.665  0.01578 * 
### ---
### Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###
### Residual standard error: 2.468 on 18 degrees of freedom
### Multiple R-squared:  0.2829,	Adjusted R-squared:  0.2431 
### F-statistic: 7.102 on 1 and 18 DF,  p-value: 0.01578


# plotten der Distanzen zur Quelle gegen log10 alpha Diversität
plot(log10(kernbach.arten.t$S)~kernbach.arten.t$dist,
     main = "Alphadiversität entlang des Kernbaches",
     xlab= "Distanz zur Quelle",
     ylab = "Anzahl der Arten")
model2<-lm(log10(kernbach.arten.t$S)~kernbach.arten.t$dist)
summary(model2)
abline(model2)
hist(model2$residuals)

# Histogramm der Residuen mit Achsen und Farbe
hist(model1$residuals,
     main = "Verteilung der Residuen",
     xlab = "Residuen",
     col = "lightblue",
     border = "white",
     breaks = 10)

# Dichtekurve ergänzen
lines(density(model1$residuals), col = "darkblue", lwd = 2)


#### Mit Berücksichtigung Fließend/Kolk sowwie der Interaktion
names(kernbach.arten.t)

model2<-lm(kernbach.arten.t$S~kernbach.arten.t$kolk*kernbach.arten.t$dist)
summary(model2)
anova(model2)



################################################################################
# I.1. Gibt es einen Zusammenhang zwischen Artenreichtum und Entfernung der 
# Probestellen von der Quelle?
# Ja.
################################################################################
################################################################################
# I.2. Unterscheiden sich stark fließenden Bereiche und die aufgestauten 
# Bereiche in Individuenreichtum und Artenreichtum?

model3<-lm(kernbach.arten.t$S~kernbach.arten.t$kolk)
summary(model3)

###Call:
###  lm(formula = kernbach.arten.t$S ~ kernbach.arten.t$kolk)
###
###Residuals:
###Min     1Q Median     3Q    Max 
###-3.6   -2.3   -0.3    1.4    6.7 
###
###Coefficients:
###  Estimate Std. Error t value Pr(>|t|)    
###(Intercept)              5.6000     0.9144   6.124 8.75e-06 ***
### kernbach.arten.t$kolkK   0.7000     1.2931   0.541    0.595    
### ---
###Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###
###Residual standard error: 2.892 on 18 degrees of freedom
###Multiple R-squared:  0.01602,	Adjusted R-squared:  -0.03865 
###F-statistic: 0.293 on 1 and 18 DF,  p-value: 0.5949


################################################################################
# I.2.a. Unterscheiden sich stark fließenden Bereiche und die aufgestauten 
# Bereiche in Artenreichtum?
# Nein - sie unterscheiden sich nicht im Artenreichtum.
################################################################################



################################################################################
################################################################################
# I.3. Ist die Artenzusammensetzung unterschiedlicher, je weiter zwei Standorte
# voneinander entfernt sind (isolation by distance)?


# - Erstellen einer Matrix der Paarweisen Unterschiede der Artenzusammensetzung
kernbach.arten.t.beta<-kernbach.arten.t
kernbach.arten.t.beta<-kernbach.arten.t.beta[,-c(44:49)]#### diese Spalten müssen auch raus weil wir für die Analyses der
#### Artenzsammensetzung nur die Matrix nutzen
kernbach.arten.t.beta<-as.data.frame(kernbach.arten.t.beta)
kernbach.arten.t.beta <- na.omit(kernbach.arten.t.beta)
str(kernbach.arten.t.beta)
### bearbeitet
daten_nur_numeric <- kernbach.arten.t.beta[, sapply(kernbach.arten.t.beta, is.numeric)]

vegdist(daten_nur_numeric, binary = TRUE, method = "jaccard")
###
beta_div_kernbach <- vegdist(daten_nur_numeric, binary = TRUE, method = "jaccard")
#### Jaccard Distanzen hier nur presence/absence
#### man kann auch andere Indices nehmen die die Abundanz der einzelnen Arten miteinberrechnen

# - Erstellen einer Matrix der paarweisen Distanzen

# 1. Nur die Artendaten extrahieren (Presence/Absence)
arten_daten <- kernbach.arten.t[, sapply(kernbach.arten.t, is.numeric)]
arten_daten <- arten_daten[, 1:38]  # Nur die Artenspalten (anpassen falls nötig)

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
     pch = 19, col = "darkblue")

# Optional: Regressionslinie einfügen
abline(lm(beta_div_vec ~ geo_dist_vec), col = "red", lwd = 2)


# Lade ggplot2, falls noch nicht geschehen
library(ggplot2)

# Erstelle das Regressionsmodell
model1 <- lm(S ~ dist, data = dist_kernbach)

# Residuen extrahieren
residuen <- model1$residuals

hist(model1$residuals)

#Residuen Copilot

# 5. Ergebnis anzeigen
print(mantel_result)

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

# Plot zur Visualisierung des Zusammenhangs
plot(df$beta ~ df$dist,
     main = "Mantel's Test: Distanz und Betadiversität",
     xlab = "Distanz zwischen zwei Aufnahmen",
     ylab = "Jaccard-Distanz")

# Lade ggplot2, falls noch nicht geschehen
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





# Regressionsmodell
model1 <- lm(beta ~ dist, data = df)

# Histogramm der Residuen
hist(model1$residuals,
     main = "Histogramm der Residuen",
     xlab = "Residuen",
     col = "lightblue",
     border = "white")


# Mantel's Test
mantel(beta_div_kernbach,dist_kernbach,  method = "pearson",
       permutations = 1000)

###Mantel statistic based on Pearson's product-moment correlation 

###Call:
###mantel(xdis = beta_div_kernbach, ydis = dist_kernbach, method = "pearson",      permutations = 1000) 
###
###Mantel statistic r: 0.3529 
###      Significance: 0.000999 
###
###Upper quantiles of permutations (null model):
###  90%   95% 97.5%   99% 
### 0.104 0.137 0.166 0.209 
### Permutation: free
### Number of permutations: 1000


################################################################################
# I.3. Ist die Artenzusammensetzung unterschiedlicher, je weiter zwei Standorte
# voneinander entfernt sind (isolation by distance)?
# Es gibt eine sigifikante positive Korrelation zwischen der Distanz und der Unähnlichkeit
# der Artenzusammensetzung der Standorte.
################################################################################
################################################################################
# 5. Ergebnis anzeigen
print(mantel_result)

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

# Visualisierung: Beta-Diversität vs. Distanz
plot(df$beta ~ df$dist,
     main = "Mantel's Test: Distanz und Betadiversität",
     xlab = "Distanz zwischen zwei Aufnahmen",
     ylab = "Jaccard-Distanz")

# Lade ggplot2
library(ggplot2)

# Modell 1: Nur Beta-Diversität
model_beta <- lm(dist ~ beta, data = df)

# Modell 2: Alpha + Beta-Diversität
model_combined <- lm(dist ~ alpha + beta, data = df)

# Modellvergleich
summary(model_beta)
summary(model_combined)
anova(model_beta, model_combined)

# Histogramm der Residuen für kombiniertes Modell
residuen <- model_combined$residuals

# Entferne NA-Werte
df_resid <- data.frame(resid = residuen)
df_resid <- na.omit(df_resid)

# Dynamische Achsenlimits – etwas Puffer für Randwerte
x_min <- floor(min(df_resid$resid)) - 1
x_max <- ceiling(max(df_resid$resid)) + 1

# Erstelle Histogramm im gewünschten Stil
ggplot(df_resid, aes(x = resid)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 2), limits = c(x_min, x_max)) +
  labs(title = "Verteilung der Residuen",
       x = "Residuen",
       y = "Häufigkeit") +
  theme_minimal(base_size = 14)

#Residuen vergleichen
# Residuen extrahieren
resid_beta <- na.omit(model_beta$residuals)
resid_combined <- na.omit(model_combined$residuals)

# Kombiniere in DataFrame
df_compare <- data.frame(
  Residuen = c(resid_beta, resid_combined),
  Modell = rep(c("Beta", "Alpha + Beta"), times = c(length(resid_beta), length(resid_combined)))
)

# Plot
ggplot(df_compare, aes(x = Residuen, fill = Modell)) +
  geom_histogram(binwidth = 1, position = "dodge", color = "black") +
  labs(title = "Vergleich der Residuen",
       x = "Residuen",
       y = "Häufigkeit") +
  theme_minimal(base_size = 14)
ggplot(df_compare, aes(x = Modell, y = Residuen, fill = Modell)) +
  geom_boxplot() +
  labs(title = "Boxplot der Residuen",
       x = "Modell",
       y = "Residuen") +
  theme_minimal(base_size = 14)

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

