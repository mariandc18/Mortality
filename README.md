# Baseline Mortality Estimation

## 📋 Descripción

El script ajusta un **Modelo Aditivo Generalizado (GAM)** con familia **Binomial Negativa** sobre los registros históricos de defunciones entre 2013 - 2020, para predecir cuántas muertes se habrían esperado en ausencia de eventos extraordinarios. La diferencia entre las muertes observadas y las esperadas constituye el exceso de mortalidad.

---

## Causas de muerte incluidas

- **Neumonía** — `J13X`, `J15`, `J154`, `J158`, `J159`, `J17`, `J170`, `J178`, `J18`, `J182`, `J188`, `J189`, `J850`, `J851`, `J440`
- **Meningitis** — `G001`, `G009`
- **Sepsis** — `A403`, `A409`, `A419`

---

## Población de estudio

Adultos mayores de 60 años, estratificados por:

| Grupo | Rango de edad |
|---|---|
| 1 | 60 – 64 años |
| 2 | 65 – 69 años |
| 3 | 70 – 74 años |
| 4 | 75 – 79 años |
| 5 | 80 – 84 años |
| 6 | 85 años o más |

---

## Modelo estadístico

Se ajusta un **GAM con familia Binomial Negativa** usando la librería `mgcv` de R:

```
observed ~ s(year, k) + s(month, bs="cc", k) + sex + age_group + offset(log(population))
```

### Componentes del modelo

- **`s(year)`** — spline suavizado sobre el año, captura la tendencia secular de la mortalidad
- **`s(month, bs="cc")`** — spline cíclico sobre el mes, captura la estacionalidad (el spline cíclico garantiza continuidad entre diciembre y enero)
- **`sex` y `age_group`** — efectos fijos por sexo y grupo de edad
- **`offset(log(population))`** — ajuste por el tamaño poblacional

---

## ⚙️ Requisitos

```r
install.packages(c("mgcv", "tidyverse", "jsonlite", "readr"))
```
