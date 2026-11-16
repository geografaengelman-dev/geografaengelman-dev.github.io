#==============================================================================
###TRABAJO DE ECONOMETRIA ESPACIAL
###Autocorrelación y Regresión espacial en R
#
###Autora: Anabella Ayelén Engelman
#==============================================================================



#==============================================================================
#PASO 1: Instalar y cargar paquetes
#==============================================================================
options(timeout = 600)
install.packages(c("spData", "spdep", "spatialreg", "tmap", 
                   "dplyr", "GWmodel", "TH.data", "lwgeom"))
install.packages("terra")
install.packages("TH.data", repos = "https://cloud.r-project.org")

#Objetos espaciales:
library(sf)
#Capa us_states:
library(spData)
#Vecindad y autocorrelación:
library(spdep)
#Modelos SAR y SEM:
library(spatialreg)
#Mapas:
library(tmap)
#Regresión geográficamente ponderada:
library(GWmodel)
#Unión de tablas y operador %>%:
library(dplyr)

setwd("C:/Users/Dell/Diplomado_26")

#==============================================================================
#PASO 2: Preparar el entorno, los datos
#==============================================================================
# 1. Datos de atributos
data(state)
edos <- as.data.frame(state.x77)

# Convertimos los nombres de las variables
names(edos) <- make.names(names(edos))

# Agregamos nombre del estado
edos$NAME <- rownames(state.x77)

# Agregamos región
edos$Region <- state.region

# Comprobamos
names(edos)
head(edos)

# 2. Geometría de los estados
us <- spData::us_states

# Comprobamos nombres de las columnas
names(us)
st_crs(us)$epsg

# 3. Verificación obligatoria: Comparación de unidades
cat("Unidades en us:", nrow(us), "\n")
cat("Unidades en edos:", nrow(edos), "\n")

cat("\nEstados de state.x77 que NO están en us_states:\n")
print(setdiff(edos$NAME, us$NAME))        #Alaska y Hawaii  

cat("\nUnidades de us_states que NO están en state.x77:\n")
print(setdiff(us$NAME, edos$NAME))        #District of Columbia


# 4. Compatibilización
#Se excluyen Alaska y Hawaii de la tabla de atributos ya que las capas
#no tienen geometría y no son contiguos a ningún estado.

edos <- edos %>%
  dplyr::filter(NAME %in% us$NAME)

# Comprobamos
cat("\nEstados con datos después del filtro:", nrow(edos), "\n")
print(edos$NAME)

# Unión de atributos y geometría
mapa <- dplyr::left_join(us, edos, by ="NAME")

cat("Unidades tras la unión:", nrow(mapa), "\n")
cat("Unidades sin datos de Murder:", sum(is.na(mapa$Murder)), "\n")
print(mapa$NAME[is.na(mapa$Murder)])     # District of Columbia


# Se excluye District of Columbia (no tiene datos en state.x77)
mapa <- mapa[!is.na(mapa$Murder), ]

# Verificación final
cat("\nUnidades sin datos después de limpiar:", sum(is.na(mapa$Murder)), "\n")
cat("Unidades finales para el análisis:", nrow(mapa), "\n")

# 5. Reproyectamos la capa
mapa <- sf::st_transform(mapa, 5070)

# 6. Modelo de Regresión
modelo <- lm(Murder ~ Illiteracy + Income + HS.Grad + Frost, data = mapa)
summary(modelo)

#==============================================================================
# PREGUNTA 1.a: Definición de la estructura de vecindad
#==============================================================================

# Centroides de cada estado (en metros, porque "mapa" ya está en EPSG:5070)
coords <- st_coordinates(st_centroid(st_geometry(mapa)))

# (i) Contigüidad: reina y torre 
nb_reina <- poly2nb(mapa)                   # queen = TRUE por defecto: frontera 
                                            # o vértice
nb_torre <- poly2nb(mapa, queen = FALSE)    # rook: solo si comparten un lado

# (ii) k vecinos más próximos, k = 4 
nb_knn <- knn2nb(knearneigh(coords, k = 4))

# (iii) Vecindad por distancia 
# Paso 1: distancia de cada estado a su vecino más cercano (con nbdists)
nb_1 <- knn2nb(knearneigh(coords, k = 1))
dist_1 <- unlist(nbdists(nb_1, coords))
summary(dist_1 / 1000)                      # en km
radio_min <- max(dist_1)                    # radio mínimo para que nadie 
                                            # quede sin vecinos
radio_min / 1000

# Paso 2: con el radio mínimo, ¿el grafo queda conectado?
nb_min <- dnearneigh(coords, 0, radio_min)
n.comp.nb(nb_min)$nc                        # número de grupos desconectados

# Paso 3: se amplía el radio hasta obtener un único grupo conectado
for (d in seq(500, 800, by = 20)) {
  cat(d, "km -> grupos:", n.comp.nb(dnearneigh(coords, 0, d * 1000))$nc, "\n")}
nb_dist <- dnearneigh(coords, 0, 580000)    # radio elegido: 580 km

# Comparación 
summary(nb_reina)
summary(nb_torre)
summary(nb_knn)
summary(nb_dist)

lista_nb <- list(Reina = nb_reina, Torre = nb_torre,
                 `KNN (k=4)` = nb_knn, `Distancia 580 km` = nb_dist)

tabla_nb <- data.frame(
  Regiones         = sapply(lista_nb, length),
  Enlaces          = sapply(lista_nb, function(n) sum(card(n))),
  Promedio         = sapply(lista_nb, function(n) round(mean(card(n)), 2)),
  Minimo           = sapply(lista_nb, function(n) min(card(n))),
  Maximo           = sapply(lista_nb, function(n) max(card(n))),
  Sin_vecinos      = sapply(lista_nb, function(n) sum(card(n) == 0)),
  Grupos_conexos   = sapply(lista_nb, function(n) n.comp.nb(n)$nc))
tabla_nb

# Estados que cambian entre reina y torre
mapa$NAME[card(nb_reina) != card(nb_torre)]

# Estados con más y con menos vecinos (reina)
mapa$NAME[card(nb_reina) == max(card(nb_reina))]
mapa$NAME[card(nb_reina) == min(card(nb_reina))]

# Gráfico: vecinos superpuestos al mapa base 
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.5, 0.5))
for (i in seq_along(lista_nb)) {
  plot(st_geometry(mapa), border = "grey60", main = names(lista_nb)[i])
  plot(lista_nb[[i]], coords, add = TRUE, col = "red")}
par(mfrow = c(1, 1))

#==============================================================================
# PREGUNTA 1.b: Matriz de pesos y autocorrelación global
#==============================================================================

# Matrices de pesos a partir de la contigüidad reina 
lw_W <- nb2listw(nb_reina, style = "W")   # estandarizada por filas
lw_B <- nb2listw(nb_reina, style = "B")   # binaria

# Ejemplo para explicar la estandarización: los pesos de Georgia
i <- which(mapa$NAME == "Georgia")
mapa$NAME[nb_reina[[i]]]                  # sus vecinos
lw_B$weights[[i]]                         # binaria: un 1 por vecino
lw_W$weights[[i]]                         # W: 1 / número de vecinos (suman 1)
sum(lw_W$weights[[i]])

# I de Moran (prueba analítica) 
moran.test(mapa$Murder, lw_W)             # Murder, pesos W
moran.test(mapa$Murder, lw_B)             # Murder, pesos B
moran.test(mapa$Illiteracy, lw_W)         # Illiteracy, pesos W

# I de Moran por permutaciones (999 simulaciones) 
set.seed(123)
mc_murder <- moran.mc(mapa$Murder, lw_W, nsim = 999)
mc_murder
set.seed(123)
mc_illit <- moran.mc(mapa$Illiteracy, lw_W, nsim = 999)
mc_illit

# Histograma de los I simulados frente al I observado
hist(mc_murder$res, breaks = 30, col = "grey85",
     main = "Moran I simulado (999 permutaciones) - Murder",
     xlab = "I de Moran", xlim = c(-0.4, 0.7))
abline(v = mc_murder$statistic, col = "red", lwd = 2)

# C de Geary 
geary.test(mapa$Murder, lw_W)
geary.test(mapa$Illiteracy, lw_W)

# Diagrama de dispersión de Moran 
moran.plot(mapa$Murder, lw_W, labels = mapa$NAME,
           xlab = "Tasa de homicidios (Murder)",
           ylab = "Rezago espacial de Murder (promedio de los vecinos)")

# Estados en cada cuadrante
z  <- mapa$Murder - mean(mapa$Murder)     # desvío respecto de la media
lz <- lag.listw(lw_W, z)                  # promedio de los desvíos de los vecinos
cuadrante <- ifelse(z > 0 & lz > 0, "Alto-Alto",
                    ifelse(z < 0 & lz < 0, "Bajo-Bajo",
                           ifelse(z > 0 & lz < 0, "Alto-Bajo", "Bajo-Alto")))
mapa$cuadrante <- cuadrante
split(mapa$NAME, cuadrante)
table(cuadrante)

#==============================================================================
# PREGUNTA 1.c: Autocorrelación local y detección de conglomerados
#==============================================================================

# Moran local (LISA) 
local <- localmoran(mapa$Murder, lw_W)
head(local)                               # Ii, E.Ii, Var.Ii, Z.Ii y valor p (columna 5)

mapa$Ii   <- local[, 1]                   # estadístico local
mapa$Ii_p <- local[, 5]                   # valor p

sum(mapa$Ii_p < 0.05)                     # cuántos estados son significativos

# Mapa del Ii y mapa de su valor p 
mapa_Ii <- tm_shape(mapa) +
  tm_fill("Ii", style = "quantile", palette = "-RdBu", midpoint = 0,
          title = "Moran local (Ii)") +
  tm_borders(col = "grey40")

mapa_p <- tm_shape(mapa) +
  tm_fill("Ii_p", breaks = c(0, 0.01, 0.05, 0.1, 1), palette = "-Greens",
          title = "Valor p") +
  tm_borders(col = "grey40")

tmap_arrange(mapa_Ii, mapa_p, ncol = 2)

# Mapa de conglomerados LISA 
# Se usa z (desvío de cada estado) y lz (promedio de los desvíos de sus
# vecinos), calculados en la 1.b.
signif <- 0.05
quadrant <- vector(mode = "numeric", length = nrow(local))
quadrant[z > 0 & lz > 0] <- 4             # alto-alto
quadrant[z < 0 & lz < 0] <- 1             # bajo-bajo
quadrant[z < 0 & lz > 0] <- 2             # bajo-alto
quadrant[z > 0 & lz < 0] <- 3             # alto-bajo
quadrant[local[, 5] > signif] <- 0        # no significativo

mapa$LISA <- factor(quadrant, levels = 0:4,
                    labels = c("No significativo", "Bajo-Bajo", "Bajo-Alto",
                               "Alto-Bajo", "Alto-Alto"))
table(mapa$LISA)

tm_shape(mapa) +
  tm_fill("LISA", palette = c("grey90", "#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"),
          title = "Conglomerados LISA") +
  tm_borders(col = "grey40")

# Estados en cada conglomerado significativo
mapa$NAME[mapa$LISA == "Alto-Alto"]
mapa$NAME[mapa$LISA == "Bajo-Bajo"]
mapa$NAME[mapa$LISA == "Alto-Bajo"]
mapa$NAME[mapa$LISA == "Bajo-Alto"]

# Getis-Ord Gi con pesos binarios 
mapa$Gi <- as.numeric(localG(mapa$Murder, lw_B))
summary(mapa$Gi)

mapa$Gi_cat <- cut(mapa$Gi,
                   breaks = c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf),
                   labels = c("Frío 99%", "Frío 95%", "Frío 90%", 
                              "No significativo",
                              "Caliente 90%", "Caliente 95%", "Caliente 99%"))
table(mapa$Gi_cat)

tm_shape(mapa) +
  tm_fill("Gi_cat", palette = c("#08519c", "#3182bd", "#9ecae1", "grey90",
                                "#fcae91", "#fb6a4a", "#cb181d"),
          title = "Getis-Ord Gi (z)") +
  tm_borders(col = "grey40")

# Puntos calientes y fríos (95%)
mapa$NAME[mapa$Gi >  1.96]                # calientes
mapa$NAME[mapa$Gi < -1.96]                # fríos

#==============================================================================
# PREGUNTA 1.d: Robustez y conclusión
#==============================================================================

# Moran de Murder con cada una de las estructuras de vecindad de la 1.a
moran.test(mapa$Murder, nb2listw(nb_reina, style = "W"))
moran.test(mapa$Murder, nb2listw(nb_torre, style = "W"))
moran.test(mapa$Murder, nb2listw(nb_knn,   style = "W"))
moran.test(mapa$Murder, nb2listw(nb_dist,  style = "W"))

# Tabla comparativa
tabla_robustez <- t(sapply(lista_nb, function(n) {
  t <- moran.test(mapa$Murder, nb2listw(n, style = "W"))
  c(I         = t$estimate[[1]],
    Esperanza = t$estimate[[2]],
    Varianza  = t$estimate[[3]],
    z         = unname(t$statistic),
    p         = t$p.value)}))
round(tabla_robustez, 4)

#==============================================================================
# GUARDAR FIGURAS Y MAPAS DE LA PREGUNTA 1
#==============================================================================

# Carpeta donde se guardarán los archivos (la carpeta raíz de trabajo)
getwd()

# --- Figura 1: grafos de vecindad (1.a) ---
png("Figura_1_grafos_vecindad.png", width = 2000, height = 1400, res = 200)
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.5, 0.5))
for (i in seq_along(lista_nb)) {
  plot(st_geometry(mapa), border = "grey60", main = names(lista_nb)[i])
  plot(lista_nb[[i]], coords, add = TRUE, col = "red")}
dev.off()

# --- Figura 2: histograma de las permutaciones de Moran (1.b) ---
png("Figura_2_moran_permutaciones.png", width = 1600, height = 1100, res = 200)
hist(mc_murder$res, breaks = 30, col = "grey85",
     main = "Moran I simulado (999 permutaciones) - Murder",
     xlab = "I de Moran", xlim = c(-0.4, 0.7))
abline(v = mc_murder$statistic, col = "red", lwd = 2)
dev.off()

# --- Figura 3: diagrama de dispersión de Moran (1.b) ---
png("Figura_3_diagrama_moran.png", width = 1600, height = 1400, res = 200)
moran.plot(mapa$Murder, lw_W, labels = mapa$NAME,
           xlab = "Tasa de homicidios (Murder)",
           ylab = "Rezago espacial de Murder (promedio de los vecinos)")
dev.off()

# --- Figura 4: Moran local, Ii y valor p (1.c) ---
tmap_save(tmap_arrange(mapa_Ii, mapa_p, ncol = 2),
          filename = "Figura_4_moran_local_Ii_p.png",
          width = 2400, height = 900, dpi = 200)

# --- Figura 5: mapa de conglomerados LISA (1.c) ---
mapa_lisa <- tm_shape(mapa) +
  tm_fill("LISA", palette = c("grey90", "#2c7bb6", "#abd9e9", "#fdae61",
                              "#d7191c"),
          title = "Conglomerados LISA") +
  tm_borders(col = "grey40")
tmap_save(mapa_lisa, filename = "Figura_5_mapa_LISA.png",
          width = 1600, height = 1000, dpi = 200)

# --- Figura 6: puntos calientes y fríos, Getis-Ord Gi (1.c) ---
mapa_gi <- tm_shape(mapa) +
  tm_fill("Gi_cat", palette = c("#08519c", "#3182bd", "#9ecae1", "grey90",
                                "#fcae91", "#fb6a4a", "#cb181d"),
          title = "Getis-Ord Gi (z)") +
  tm_borders(col = "grey40")
tmap_save(mapa_gi, filename = "Figura_6_getis_ord.png",
          width = 1600, height = 1000, dpi = 200)

# Comprobación: lista de figuras guardadas
list.files(pattern = "^Figura_.*\\.png$")

#==============================================================================
# PREGUNTA 2.a: El modelo MCO y sus residuos
#==============================================================================

# Resumen del modelo de referencia 
summary(modelo)

# Los cuatro gráficos de diagnóstico 
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(modelo, labels.id = mapa$NAME)
par(mfrow = c(1, 1))

# Guardar la figura
png("Figura_7_diagnostico_MCO.png", width = 2000, height = 1800, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(modelo, labels.id = mapa$NAME)
dev.off()

# Residuos y residuos estandarizados 
mapa$res  <- residuals(modelo)            # residuo = valor real - valor predicho
mapa$rstd <- rstandard(modelo)            # residuo estandarizado

summary(mapa$rstd)

# Estados con los residuos más grandes
head(st_drop_geometry(mapa)[order(-mapa$rstd), c("NAME", "Murder", "rstd")], 4)   
# subestimados por el modelo
head(st_drop_geometry(mapa)[order(mapa$rstd),  c("NAME", "Murder", "rstd")], 4)   
# sobreestimados por el modelo

# Mapa de residuos con paleta divergente 
mapa_res <- tm_shape(mapa) +
  tm_fill("rstd", palette = "-RdBu", midpoint = 0, style = "fixed",
          breaks = c(-3, -2, -1, -0.5, 0.5, 1, 2, 3.5),
          title = "Residuo estandarizado") +
  tm_borders(col = "grey40")
mapa_res

tmap_save(mapa_res, filename = "Figura_8_mapa_residuos_MCO.png",
          width = 1600, height = 1000, dpi = 200)

# I de Moran sobre los residuos del MCO 
moran_res_mco <- lm.morantest(modelo, lw_W)
moran_res_mco

#==============================================================================
# PREGUNTA 2.b: Diagnóstico de especificación (tests LM)
#==============================================================================

# Los cinco tests de multiplicadores de Lagrange
tests_lm <- lm.LMtests(modelo, lw_W, test = "all")
tests_lm

# Tabla resumen: estadístico, grados de libertad y valor p
tabla_lm <- data.frame(
  Estadistico = sapply(tests_lm, function(x) unname(x$statistic)),
  gl          = sapply(tests_lm, function(x) unname(x$parameter)),
  p           = sapply(tests_lm, function(x) x$p.value))
round(tabla_lm, 4)

#==============================================================================
# PREGUNTA 2.c: Estimación, comparación y validación de los modelos globales
#==============================================================================

# Modelo de rezago espacial (SAR) 
m_lag <- lagsarlm(Murder ~ Illiteracy + Income + HS.Grad + Frost,
                  data = mapa, listw = lw_W)
summary(m_lag)

# Modelo de error espacial (SEM) 
m_err <- errorsarlm(Murder ~ Illiteracy + Income + HS.Grad + Frost,
                    data = mapa, listw = lw_W)
summary(m_err)

# Test de razón de verosimilitud contra el MCO 
LR.Sarlm(m_lag, modelo)                   # SAR frente a MCO
LR.Sarlm(m_err, modelo)                   # SEM frente a MCO

# Tabla comparativa de los tres modelos 
s_lag <- summary(m_lag)
s_err <- summary(m_err)
lr_lag <- LR.Sarlm(m_lag, modelo)
lr_err <- LR.Sarlm(m_err, modelo)

tabla_modelos <- data.frame(
  MCO = c(coef(modelo)["Illiteracy"], coef(modelo)["Income"],
          coef(modelo)["HS.Grad"], coef(modelo)["Frost"],
          NA, NA, logLik(modelo), AIC(modelo), NA, NA),
  SAR = c(coef(m_lag)["Illiteracy"], coef(m_lag)["Income"],
          coef(m_lag)["HS.Grad"], coef(m_lag)["Frost"],
          m_lag$rho, s_lag$Wald1$p.value, logLik(m_lag), AIC(m_lag),
          lr_lag$statistic, lr_lag$p.value),
  SEM = c(coef(m_err)["Illiteracy"], coef(m_err)["Income"],
          coef(m_err)["HS.Grad"], coef(m_err)["Frost"],
          m_err$lambda, s_err$Wald1$p.value, logLik(m_err), AIC(m_err),
          lr_err$statistic, lr_err$p.value),
  row.names = c("Illiteracy", "Income", "HS.Grad", "Frost",
                "rho / lambda", "p (rho / lambda)", "Log-verosimilitud",
                "AIC", "LR vs MCO", "p (LR)"))
round(tabla_modelos, 4)

# Efectos directos, indirectos y totales (modelo SAR)
set.seed(123)
impactos <- impacts(m_lag, listw = lw_W, R = 1000)
summary(impactos, zstats = TRUE, short = TRUE)

# Validación: Moran sobre los residuos del SAR
moran_res_sar <- moran.test(residuals(m_lag), lw_W)
moran_res_sar

# Comparación con los residuos del MCO (2.a)
data.frame(Modelo  = c("MCO", "SAR"),
           I_Moran = c(moran_res_mco$estimate[[1]], moran_res_sar$estimate[[1]]),
           p       = c(moran_res_mco$p.value,       moran_res_sar$p.value))

#==============================================================================
# PREGUNTA 2.d: Regresión geográficamente ponderada (GWR)
#==============================================================================

mapa_sp <- as(mapa, "Spatial")

# Selección del ancho de banda óptimo 

# kernel bicuadrado y ancho de banda adaptativo (número de vecinos), criterio AICc
bw <- bw.gwr(Murder ~ Illiteracy + Income + HS.Grad + Frost,
             data = mapa_sp, approach = "AICc",
             kernel = "bisquare", adaptive = TRUE)
bw                                  # número de vecinos en cada regresión local

# Estimación de la GWR 
gwr_modelo <- gwr.basic(Murder ~ Illiteracy + Income + HS.Grad + Frost,
                        data = mapa_sp, bw = bw,
                        kernel = "bisquare", adaptive = TRUE)
gwr_modelo

#  Resultados locales 
mapa$b_Illit  <- gwr_modelo$SDF$Illiteracy  # coeficiente local del analfabetismo
mapa$t_Illit  <- gwr_modelo$SDF$Illiteracy_TV # su valor t local
mapa$R2_local <- gwr_modelo$SDF$Local_R2    # R2 local

summary(mapa$b_Illit)
summary(mapa$R2_local)

# Estados donde el efecto del analfabetismo es mayor y menor
head(st_drop_geometry(mapa)[order(-mapa$b_Illit), 
                            c("NAME", "b_Illit", "t_Illit", "R2_local")], 6)
head(st_drop_geometry(mapa)[order(mapa$b_Illit),  
                            c("NAME", "b_Illit", "t_Illit", "R2_local")], 6)

# Estados donde el efecto es significativo (|t| > 1,96)
mapa$NAME[abs(mapa$t_Illit) > 1.96]

# Comparación de ajuste GWR vs MCO
c(R2_ajustado_MCO = summary(modelo)$adj.r.squared,
  R2_ajustado_GWR = gwr_modelo$GW.diagnostic$gwR2.adj,
  AICc_GWR        = gwr_modelo$GW.diagnostic$AICc)

# Mapas: coeficiente de Illiteracy y R2 local, en una sola figura 
mapa_coef <- tm_shape(mapa) +
  tm_fill("b_Illit", palette = "-RdBu", midpoint = 0, style = "fixed",
          breaks = c(-3.1, -2, -1, 0, 1, 2, 3, 4, 5, 5.6),
          title = "Coef. local Illiteracy") +
  tm_borders(col = "grey40")

mapa_r2 <- tm_shape(mapa) +
  tm_fill("R2_local", palette = "YlGn", style = "quantile",
          title = "R2 local") +
  tm_borders(col = "grey40")

tmap_arrange(mapa_coef, mapa_r2, ncol = 2)

tmap_save(tmap_arrange(mapa_coef, mapa_r2, ncol = 2),
          filename = "Figura_9_GWR_coef_R2.png",
          width = 2400, height = 900, dpi = 200)

#==============================================================================
# GUARDAMOS LOS PRODUCTOS DE LA PREGUNTA 2
#==============================================================================

getwd()                                   

#------------------------------------------------------------------------------
# 1. FIGURAS
#------------------------------------------------------------------------------

# Figura 7: diagnóstico del MCO (2.a)
png("Figura_7_diagnostico_MCO.png", width = 2000, height = 1800, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(modelo, labels.id = mapa$NAME)
dev.off()

# Figura 8: mapa de residuos del MCO (2.a)
tmap_save(mapa_res, filename = "Figura_8_mapa_residuos_MCO.png",
          width = 1600, height = 1000, dpi = 200)

# Figura 9: GWR, coeficiente local de Illiteracy y R2 local (2.d)
tmap_save(tmap_arrange(mapa_coef, mapa_r2, ncol = 2),
          filename = "Figura_9_GWR_coef_R2.png",
          width = 2400, height = 900, dpi = 200)

#------------------------------------------------------------------------------
# 2. TABLAS 
#------------------------------------------------------------------------------

# Tabla 4: tests de multiplicadores de Lagrange (2.b)
write.csv2(round(tabla_lm, 4), "Tabla_4_tests_LM.csv")

# Tabla 5: comparación de los modelos MCO, SAR y SEM (2.c)
write.csv2(round(tabla_modelos, 4), "Tabla_5_comparacion_modelos.csv")

# Tabla 6: efectos directos, indirectos y totales del SAR (2.c)
res_imp <- summary(impactos, zstats = TRUE, short = TRUE)
tabla_impactos <- data.frame(
  Directo     = impactos$res$direct,
  p_Directo   = res_imp$pzmat[, "Direct"],
  Indirecto   = impactos$res$indirect,
  p_Indirecto = res_imp$pzmat[, "Indirect"],
  Total       = impactos$res$total,
  p_Total     = res_imp$pzmat[, "Total"],
  row.names   = attr(impactos, "bnames"))
round(tabla_impactos, 4)
write.csv2(round(tabla_impactos, 4), "Tabla_6_efectos_SAR.csv")

# Tabla 7: Moran de los residuos, MCO frente a SAR (2.c)
tabla_moran_res <- data.frame(
  Modelo    = c("MCO", "SAR"),
  I_Moran   = c(moran_res_mco$estimate[[1]], moran_res_sar$estimate[[1]]),
  Esperanza = c(moran_res_mco$estimate[[2]], moran_res_sar$estimate[[2]]),
  p         = c(moran_res_mco$p.value,       moran_res_sar$p.value))
write.csv2(tabla_moran_res, "Tabla_7_moran_residuos.csv", row.names = FALSE)

# Tabla 8: resultados locales de la GWR por estado (2.d)
tabla_gwr <- st_drop_geometry(mapa)[, c("NAME", "b_Illit", "t_Illit", "R2_local")]
tabla_gwr <- tabla_gwr[order(-tabla_gwr$b_Illit), ]
write.csv2(tabla_gwr, "Tabla_8_GWR_por_estado.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# 3. SALIDAS COMPLETAS PARA EL ANEXO
#------------------------------------------------------------------------------
sink("Anexo_Pregunta_2_salidas.txt")

cat("=== 2.a Modelo MCO ===\n");                     print(summary(modelo))
cat("\n=== 2.a Moran sobre residuos del MCO ===\n");  print(moran_res_mco)
cat("\n=== 2.b Tests LM ===\n");                      print(tests_lm)
cat("\n=== 2.c Modelo SAR ===\n");                    print(summary(m_lag))
cat("\n=== 2.c Modelo SEM ===\n");                    print(summary(m_err))
cat("\n=== 2.c LR SAR vs MCO ===\n");                 print(LR.Sarlm(m_lag, 
                                                                     modelo))
cat("\n=== 2.c LR SEM vs MCO ===\n");                 print(LR.Sarlm(m_err, 
                                                                     modelo))
cat("\n=== 2.c Efectos del SAR ===\n");               print(res_imp)
cat("\n=== 2.c Moran sobre residuos del SAR ===\n");  print(moran_res_sar)
cat("\n=== 2.d Ancho de banda GWR ===\n");            print(bw)
cat("\n=== 2.d GWR ===\n");                           print(gwr_modelo)

sink()

#------------------------------------------------------------------------------
# 4. COMPROBACIÓN: lista de archivos guardados
#------------------------------------------------------------------------------
list.files(pattern = "^(Figura_[789]|Tabla_[4-8]|Anexo_Pregunta_2)")

#==============================================================================
# INFORMACIÓN DE LA SESIÓN
#==============================================================================
sessionInfo()
