# 1. LIBRERÍAS ---------------------------------------------------------------
library(ggplot2)
library(ggfortify)
library(aTSA)
library(nortest)
library(moments)
library(fitdistrplus)
library(MASS)
library(tseries)
library(tsoutliers)
library(forecast)
library(expsmooth)
library(fma)
library(FinTS)
library(rugarch)
library(scales)
library(ggpubr)
#2. DATOS ----------------------------------------------------------------------
url="https://raw.githubusercontent.com/Zairpv/Time-series-modeling---Climate/refs/heads/main/Data/Precipitacion.csv"
Precipit <- read.csv(url, header = TRUE)
View(Precipit)
# 3. SERIE DE TIEMPO -----------------------------------------------------------
precipitacion <- as.numeric(t(as.matrix(Precipit[, -1])))
serie_precipitacion <- ts(precipitacion, start = c(1960, 1), frequency = 12)
plot(serie_precipitacion,
     main = "Serie Temporal: Precipitación del monte Tláloc",
     ylab = "Precipitación (mm)",
     xlab = "Años")

# Extrae ventana de interés: 2002 a 2017
Periodo_2002_2017 <- window(serie_precipitacion, start = c(2002, 1), end = c(2017, 12))

# 4. ATRIBUTOS DE LA SERIE ------------------------------------------------------

# Visualización básica de la serie recortada
plot(Periodo_2002_2017,
     main = "Precipitación del monte Tláloc (2002–2017)",
     ylab = "Precipitación (mm)",
     xlab = "Años")

# a) Descomposición clásica: tendencia, estacionalidad y componente aleatorio
descomposición <- decompose(Periodo_2002_2017)
plot(descomposición)
# Notas:
# - Trend: evolución a largo plazo (↑ o ↓ en precipitación)
# - Seasonal: patrón promedio mensual
# - Random: variación no explicada por los dos anteriores

# b) Funciones de autocorrelación (ACF) y autocorrelación parcial (PACF)
par(mfrow = c(2, 1), mar = c(5, 4, 3, 1), cex = 0.6)
acf(Periodo_2002_2017,
    main = "Autocorrelación (ACF)")
pacf(Periodo_2002_2017,
     main = "Autocorrelación Parcial (PACF)")
dev.off()


# 5. ATRIBUTOS DE LA SERIE TRANSFORMADA -----------------------------------------

# Aplicación de transformación logarítmica para estabilizar la varianza
Periodo_2002_2017_log <- log(Periodo_2002_2017 + 1)

# a) Histogramas comparativos: sin transformar vs log-transformada
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1), cex = 0.6)
hist(Periodo_2002_2017,
     breaks = 20,
     col = "lightgrey",
     freq = FALSE,
     main = "Serie sin transformar",
     xlab = "Precipitación (mm)")
lines(density(Periodo_2002_2017), col = "red", lwd = 2)

hist(Periodo_2002_2017_log,
     breaks = 20,
     col = "lightgrey",
     freq = FALSE,
     main = "Serie transformada: log(1 + x)",
     xlab = "log(1 + precip)")
lines(density(Periodo_2002_2017_log), col = "red", lwd = 2)
dev.off()

# b) Estadísticos descriptivos y pruebas de normalidad
(media_log <- mean(Periodo_2002_2017_log))
(mediana_log <- median(Periodo_2002_2017_log))
(asimetria_log <- skewness(Periodo_2002_2017_log))
(curtosis_log <- kurtosis(Periodo_2002_2017_log))
(jb_test_log <- jarque.bera.test(Periodo_2002_2017_log))

# c) Descomposición de la serie transformada
descomposición_log <- decompose(Periodo_2002_2017_log)

# Nota: La transformación no altera la estacionalidad, pero mejora la estabilidad de los residuales
plot(descomposición_log)

# d) ACF y PACF de la serie transformada
par(mfrow = c(2, 1), mar = c(5, 4, 3, 1), cex = 0.6)
acf(Periodo_2002_2017_log,
    main = "ACF: log(1 + precipitación)")
pacf(Periodo_2002_2017_log,
     main = "PACF: log(1 + precipitación)")
dev.off()

# Evaluación adicional de autocorrelación
# El gráfico generado por identify() muestra los p-values del test de Ljung-Box para diferentes lags
identify(Periodo_2002_2017_log, stat.test = TRUE)
# Nota:
# - p-values

#6. MODELOS ARIMA ---------------------------------------------------------------

## 6.1. Evaluación de ordenes del modelo
# Evaluación automática de órdenes para comparar modelos con y sin transformación
auto.arima(Periodo_2002_2017_log)    # → (1,0,0)(1,1,0)[12]


## 6.2. Modelo ARIMA simple (log-transformado)
ari.m1 <- Arima(Periodo_2002_2017_log,
                order = c(1, 0, 0),
                seasonal = list(order = c(1, 1, 0), period = 12),
                method = "ML")

summary(ari.m1)
#Notas
#AR(1): efecto de corto plazo leve.
#SAR(1): captura clara de estacionalidad.


## 6.3 Diagnóstico de residuales (ARIMA simple)
checkresiduals(ari.m1)
Box.test(residuals(ari.m1), lag = 24, fitdf = 2, type = "Ljung")  # Autocorrelación
ArchTest(residuals(ari.m1), lags = 12)                            # Heterocedasticidad
#La serie no es completamente explicada, lo que motiva el análisis de intervenciones.


## 6.4. Identificación de outliers estructurales
resid_arima1 <- residuals(ari.m1)
pars_arima1 <- coefs2poly(ari.m1)

valores_atipicos_M1 <- locate.outliers(resid_arima1, pars_arima1)
valores_atipicos_M1

# Se detectan los siguientes eventos:
time(Periodo_2002_2017_log)[c(35, 59, 62, 99, 122)]
#Tipos:
#AO = outlier aditivo
#TC = cambio temporal (efecto gradual)

# Matriz de intervención
intervenciones_M1 <- matrix(0, nrow = length(Periodo_2002_2017_log), ncol = 5)
colnames(intervenciones_M1) <- c("AO_62", "TC_35", "TC_59", "TC_99", "TC_122")

intervenciones_M1[62, 1] <- 1
intervenciones_M1[35, 2] <- 1
intervenciones_M1[59, 3] <- 1
intervenciones_M1[99, 4] <- 1
intervenciones_M1[122, 5] <- 1

#Visualización de los outliers en la serie
tiempo <- seq(as.Date("2002-01-01"), by = "month", length.out = 192)
outliers_indices <- c(35, 59, 62, 99, 122)

plot(tiempo, Periodo_2002_2017_log, type = "l", col = "black", lwd = 2,
     main = "Precipitaciones log-transformadas (2002–2017)",
     xlab = "Tiempo", ylab = "Log(Precipitación)")
points(tiempo[outliers_indices], Periodo_2002_2017_log[outliers_indices],
       col = "red", pch = 19, cex = 1.3)
legend("topright", legend = "Outliers", col = "red", pch = 19)


## 6.5 Modelo ARIMA + intervención
arim1_interv <- Arima(Periodo_2002_2017_log,
                      order = c(1, 0, 0),
                      seasonal = list(order = c(1, 1, 0), period = 12),
                      xreg = intervenciones_M1,
                      method = "ML")

summary(arim1_interv)

#Evaluación de coeficientes
coefs <- coef(arim1_interv)
se <- sqrt(diag(vcov(arim1_interv)))
t_values <- coefs / se
p_values <- 2 * (1 - pnorm(abs(t_values)))

resultados <- data.frame(
  Estimación = coefs,
  Error_Estandar = se,
  t_value = t_values,
  p_value = p_values
)
round(resultados, 4)
#TC_99, AO_62 y TC_122 son significativos (p < 0.01)

#Diagnóstico del modelo con intervención
checkresiduals(arim1_interv)
Box.test(residuals(arim1_interv), lag = 24, fitdf = 7, type = "Ljung")
ArchTest(residuals(arim1_interv), lags = 12)


## 6.6. Comparación de modelos (con vs sin intervención)
AIC(arim1_interv, ari.m1)
BIC(arim1_interv, ari.m1)

#Visualización: log(1 + precip)
ajustado <- fitted(arim1_interv)

df_log <- data.frame(
  Fecha = tiempo,
  Observado = as.numeric(Periodo_2002_2017_log),
  Ajustado = as.numeric(ajustado)
)

ggplot(df_log, aes(x = Fecha)) +
  geom_line(aes(y = Observado, color = "Observado"), size = 1) +
  geom_line(aes(y = Ajustado, color = "Ajustado"), size = 1) +
  scale_color_manual(values = c("Observado" = "black", "Ajustado" = "blue")) +
  labs(title = "Precipitación log-transformada vs modelo con intervención",
       y = "Log(Precipitación + 1)", x = "Tiempo", color = "") +
  theme_minimal() +
  theme(legend.position = "top")


# Visualización: escala original
Serie_original <- exp(Periodo_2002_2017_log) - 1
ajustado_original <- exp(ajustado) - 1

df_original <- data.frame(
  Fecha = tiempo,
  Observado = Serie_original,
  Ajustado = ajustado_original
)

ggplot(df_original, aes(x = Fecha)) +
  geom_line(aes(y = Observado, color = "Observado"), size = 1) +
  geom_line(aes(y = Ajustado, color = "Ajustado"), size = 0.8, alpha = 0.7) +
  geom_point(data = df_original[outliers_indices, ],
             aes(x = Fecha, y = Observado), color = "red", size = 2) +
  scale_color_manual(values = c("Observado" = "black", "Ajustado" = "blue")) +
  labs(title = "Precipitación mensual (escala original)",
       y = "Precipitación (mm)", x = "Tiempo", color = "") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2017-12-01"))) +
  theme(legend.position = "top")


#7. PREDICCIONES ---------------------------------------------------------------
##7.1. PARA MODELO ARIMA SIMPLE
# Pronóstico de 12 meses con el modelo sin intervención
pred_simple <- forecast::forecast(ari.m1, h = 12)


# Retransformar a escala original
# Aplicar exp() - 1 a la predicción y a los intervalos
pred_simple_mean <- exp(pred_simple$mean) - 1
pred_simple_lower <- exp(pred_simple$lower[,2]) - 1  # intervalo 95%
pred_simple_upper <- exp(pred_simple$upper[,2]) - 1

# Fechas del pronóstico
tiempo_pred_simple <- seq(max(tiempo) + 1, by = "month", length.out = 12)

# Data frame combinado
df_pred_simple <- data.frame(
  Fecha = c(tiempo, tiempo_pred_simple),
  Observado = c(Serie_original, rep(NA, 12)),
  Ajustado = c(ajustado_original, pred_simple_mean),
  Lower = c(rep(NA, length(tiempo)), pred_simple_lower),
  Upper = c(rep(NA, length(tiempo)), pred_simple_upper)
)

PREDICCIONMODELOARIMA_SIMPLE=ggplot(df_pred_simple, aes(x = Fecha)) +
  geom_line(aes(y = Observado), color = "black", size = 1) +
  geom_line(aes(y = Ajustado), color = "blue", size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  scale_y_continuous(limits = c(0,600),breaks = c(0,100,200,300,400,500,600))+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2018-12-01"))) +
  labs(title = "Predicción: ARIMA (1,0,0) (1,1,0) [12]",
       subtitle = "Horizonte de 12 meses con intervalo de confianza al 95%",
       x = "Tiempo", y = "Precipitación (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))
PREDICCIONMODELOARIMA_SIMPLE

##7.2. PARA MODELO ARIMA + INTERVENCIÓN -----------------------------------------
# Crear matriz de 12 filas (1 año) y 5 columnas (una por intervención usada)
xreg_futuro <- matrix(0, nrow = 12, ncol = 5)
colnames(xreg_futuro) <- colnames(intervenciones_M1)
# Predecir 12 pasos adelante
prediccion <- forecast(arim1_interv, xreg = xreg_futuro, h = 12)

# Ver resumen
summary(prediccion)

# Convertir predicción y límites a escala original
pred_original <- exp(prediccion$mean) - 1
lower_original <- exp(prediccion$lower[,2]) - 1  # intervalo 95%
upper_original <- exp(prediccion$upper[,2]) - 1

# Crear fechas para predicción
tiempo_pred <- seq(max(tiempo) + 1, by = "month", length.out = 12)

# Data frame con observaciones y predicción
df_pred <- data.frame(
  Fecha = c(tiempo, tiempo_pred),
  Observado = c(Serie_original, rep(NA, 12)),
  Ajustado = c(ajustado_original, pred_original),
  Lower = c(rep(NA, length(tiempo)), lower_original),
  Upper = c(rep(NA, length(tiempo)), upper_original)
)

# Graficar
PREDICCIONMODELOARIMA_INTERVENCION=ggplot(df_pred, aes(x = Fecha)) +
  geom_line(aes(y = Observado), color = "black", size = 1) +
  geom_line(aes(y = Ajustado), color = "blue", size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  geom_point(data = df_pred[outliers_indices, ],
             aes(x = Fecha, y = Observado), color = "red", size = 2) +
  labs(title = "Predicción: ARIMA (1,0,0) (1,1,0) [12] + Intervención",
       subtitle = "Horizonte de 12 meses con intervalo de confianza al 95%",
       x = "Tiempo", y = "Precipitación (mm)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2018-12-01"))) +
  scale_y_continuous(limits = c(0,600),breaks = c(0,100,200,300,400,500,600))+
  theme_minimal()
PREDICCIONMODELOARIMA_INTERVENCION

#tabla con los valores pronosticados:
pred_tabla <- data.frame(
  Fecha = tiempo_pred,
  Prediccion = round(pred_original, 1),
  LimiteInferior = round(lower_original, 1),
  LimiteSuperior = round(upper_original, 1)
)

print(pred_tabla)


##7.3. AMBOS MODELOS EN ESCALA ORIGINAL -----------------------------------------------------------
ggarrange(PREDICCIONMODELOARIMA_SIMPLE,
          PREDICCIONMODELOARIMA_INTERVENCION,
          ncol = 1,
          nrow = 2)

