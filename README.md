# Time-series-modeling---Climate
Creation of ARIMA model with interventions for precipitation data

# 1. Introducción

Este informe documenta el análisis de series de tiempo aplicado a datos mensuales de precipitación en Texcoco (2002–2017), con énfasis en la modelación ARIMA y la identificación de eventos de intervención para mejorar la capacidad predictiva.

# 2. Librerías y configuración

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
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
```

# 3. Carga de datos

```{r datos}
url="https://raw.githubusercontent.com/Zairpv/Time-series-modeling---Climate/main/Data/Precipitacion.csv"
Precipit <- read.csv(url, header = TRUE)
precipitacion <- as.numeric(t(as.matrix(Precipit[, -1])))
serie_precipitacion <- ts(precipitacion, start = c(1960, 1), frequency = 12)
Periodo_2002_2017 <- window(serie_precipitacion, start = c(2002, 1), end = c(2017, 12))
```

# 4. Exploración de la serie de tiempo

```{r exploracion-serie}
plot(serie_precipitacion, main = "Serie Temporal: Precipitación del monte Tláloc",
     ylab = "Precipitación (mm)", xlab = "Años")

plot(Periodo_2002_2017,
     main = "Precipitación del monte Tláloc (2002–2017)",
     ylab = "Precipitación (mm)",
     xlab = "Años")

descomposición <- decompose(Periodo_2002_2017)
plot(descomposición)

par(mfrow = c(2, 1), mar = c(5, 4, 3, 1), cex = 0.6)
acf(Periodo_2002_2017, main = "Autocorrelación (ACF)")
pacf(Periodo_2002_2017, main = "Autocorrelación Parcial (PACF)")
dev.off()
```

# 5. Transformación logarítmica y análisis estadístico

```{r transformacion}
Periodo_2002_2017_log <- log(Periodo_2002_2017 + 1)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 1), cex = 0.6)
hist(Periodo_2002_2017, breaks = 20, col = "lightgrey", freq = FALSE,
     main = "Serie sin transformar", xlab = "Precipitación (mm)")
lines(density(Periodo_2002_2017), col = "red", lwd = 2)

hist(Periodo_2002_2017_log, breaks = 20, col = "lightgrey", freq = FALSE,
     main = "Serie transformada: log(1 + x)", xlab = "log(1 + precip)")
lines(density(Periodo_2002_2017_log), col = "red", lwd = 2)
dev.off()

(media_log <- mean(Periodo_2002_2017_log))
(mediana_log <- median(Periodo_2002_2017_log))
(asimetria_log <- skewness(Periodo_2002_2017_log))
(curtosis_log <- kurtosis(Periodo_2002_2017_log))
(jb_test_log <- jarque.bera.test(Periodo_2002_2017_log))
```

```{r acf-pacf-log}
descomposición_log <- decompose(Periodo_2002_2017_log)
plot(descomposición_log)

par(mfrow = c(2, 1), mar = c(5, 4, 3, 1), cex = 0.6)
acf(Periodo_2002_2017_log, main = "ACF: log(1 + precipitación)")
pacf(Periodo_2002_2017_log, main = "PACF: log(1 + precipitación)")
dev.off()
```

# 6. Modelado ARIMA

```{r arima-modelo}
auto.arima(Periodo_2002_2017_log)
ari.m1 <- Arima(Periodo_2002_2017_log,
                order = c(1, 0, 0),
                seasonal = list(order = c(1, 1, 0), period = 12),
                method = "ML")
summary(ari.m1)
checkresiduals(ari.m1)
Box.test(residuals(ari.m1), lag = 24, fitdf = 2, type = "Ljung")
ArchTest(residuals(ari.m1), lags = 12)
```

# 7. Identificación de intervenciones y modelo ajustado

```{r intervencion}
resid_arima1 <- residuals(ari.m1)
pars_arima1 <- coefs2poly(ari.m1)
valores_atipicos_M1 <- locate.outliers(resid_arima1, pars_arima1)
valores_atipicos_M1
time(Periodo_2002_2017_log)[c(35, 59, 62, 99, 122)]

intervenciones_M1 <- matrix(0, nrow = length(Periodo_2002_2017_log), ncol = 5)
colnames(intervenciones_M1) <- c("AO_62", "TC_35", "TC_59", "TC_99", "TC_122")
intervenciones_M1[62, 1] <- 1
intervenciones_M1[35, 2] <- 1
intervenciones_M1[59, 3] <- 1
intervenciones_M1[99, 4] <- 1
intervenciones_M1[122, 5] <- 1

tiempo <- seq(as.Date("2002-01-01"), by = "month", length.out = 192)
outliers_indices <- c(35, 59, 62, 99, 122)

plot(tiempo, Periodo_2002_2017_log, type = "l", col = "black", lwd = 2,
     main = "Precipitaciones log-transformadas (2002–2017)",
     xlab = "Tiempo", ylab = "Log(Precipitación)")
points(tiempo[outliers_indices], Periodo_2002_2017_log[outliers_indices],
       col = "red", pch = 19, cex = 1.3)
legend("topright", legend = "Outliers", col = "red", pch = 19)

arim1_interv <- Arima(Periodo_2002_2017_log,
                      order = c(1, 0, 0),
                      seasonal = list(order = c(1, 1, 0), period = 12),
                      xreg = intervenciones_M1,
                      method = "ML")
summary(arim1_interv)
```

```{r comparacion-ajuste}
ajustado <- fitted(arim1_interv)
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
  scale_color_manual(values = c("Observado" = "black", "Ajustado" = "blue")) +
  labs(title = "Precipitación mensual (escala original)",
       y = "Precipitación (mm)", x = "Tiempo", color = "") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2017-12-01"))) +
  theme(legend.position = "top")
```

# 8. Predicciones (12 meses)

```{r predicciones}
pred_simple <- forecast::forecast(ari.m1, h = 12)
pred_simple_mean <- exp(pred_simple$mean) - 1
pred_simple_lower <- exp(pred_simple$lower[,2]) - 1
pred_simple_upper <- exp(pred_simple$upper[,2]) - 1
tiempo_pred_simple <- seq(max(tiempo) + 1, by = "month", length.out = 12)

df_pred_simple <- data.frame(
  Fecha = c(tiempo, tiempo_pred_simple),
  Observado = c(Serie_original, rep(NA, 12)),
  Ajustado = c(ajustado_original, pred_simple_mean),
  Lower = c(rep(NA, length(tiempo)), pred_simple_lower),
  Upper = c(rep(NA, length(tiempo)), pred_simple_upper)
)

PREDICCIONMODELOARIMA_SIMPLE <- ggplot(df_pred_simple, aes(x = Fecha)) +
  geom_line(aes(y = Observado), color = "black", size = 1) +
  geom_line(aes(y = Ajustado), color = "blue", size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  scale_y_continuous(limits = c(0,600),breaks = seq(0,600,100)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2018-12-01"))) +
  labs(title = "Predicción: ARIMA (1,0,0)(1,1,0)[12]",
       subtitle = "12 meses con intervalo de confianza 95%",
       x = "Tiempo", y = "Precipitación (mm)") +
  theme_minimal()

PREDICCIONMODELOARIMA_SIMPLE
```

```{r prediccion-intervencion}
xreg_futuro <- matrix(0, nrow = 12, ncol = 5)
colnames(xreg_futuro) <- colnames(intervenciones_M1)
prediccion <- forecast(arim1_interv, xreg = xreg_futuro, h = 12)

pred_original <- exp(prediccion$mean) - 1
lower_original <- exp(prediccion$lower[,2]) - 1
upper_original <- exp(prediccion$upper[,2]) - 1
tiempo_pred <- seq(max(tiempo) + 1, by = "month", length.out = 12)

df_pred <- data.frame(
  Fecha = c(tiempo, tiempo_pred),
  Observado = c(Serie_original, rep(NA, 12)),
  Ajustado = c(ajustado_original, pred_original),
  Lower = c(rep(NA, length(tiempo)), lower_original),
  Upper = c(rep(NA, length(tiempo)), upper_original)
)

PREDICCIONMODELOARIMA_INTERVENCION <- ggplot(df_pred, aes(x = Fecha)) +
  geom_line(aes(y = Observado), color = "black", size = 1) +
  geom_line(aes(y = Ajustado), color = "blue", size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  geom_point(data = df_pred[outliers_indices, ],
             aes(x = Fecha, y = Observado), color = "red", size = 2) +
  labs(title = "Predicción: ARIMA (1,0,0)(1,1,0)[12] + Intervención",
       subtitle = "12 meses con intervalo de confianza 95%",
       x = "Tiempo", y = "Precipitación (mm)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = as.Date(c("2002-01-01", "2018-12-01"))) +
  scale_y_continuous(limits = c(0,600),breaks = seq(0,600,100)) +
  theme_minimal()

PREDICCIONMODELOARIMA_INTERVENCION
```

```{r comparacion-final}
ggarrange(PREDICCIONMODELOARIMA_SIMPLE,
          PREDICCIONMODELOARIMA_INTERVENCION,
          ncol = 1, nrow = 2)
```

