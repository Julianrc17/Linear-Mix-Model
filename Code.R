library(readxl)
library(nlme)
library(lme4)
library(lmerTest)
library(Matrix)
library(interactions)
#library(ggiraph)
#library(ggiraphExtra)
library(plyr)
library(jtools)
library(ggplot2)
library(ggeffects)
library(dplyr)
library(emmeans)
library(performance)
library(sjPlot)
library(knitr)
library(rmarkdown)
library(kableExtra)

df <- read.table("FACEHBI_BD_long_Araclon.ver2.0.txt", sep="\t", header=TRUE)
dfl <- read.table("FACEHBI_BD_long_Araclon.ver2.0.long.txt", sep="\t", header=TRUE)
colnames(dfl)
dim(dfl)
sum(is.na(dfl$Time))
sum(is.na(dfl$facehbi_id_fac))
class(dfl$apoe4_al)
class(dfl$Time)
class(df$Plasma_AB42_40_V0_dic)

lista_variables <- c('fbb_pet_centilod', 'aHipocampal_vol', 'eTIV', 'aVentricular_vol', 'Cortical_vol', 
                     'mmse_fac', 'fname_total_fac', 'fe_skt_time_np', 'rbans_copy_fac', 'rbans_delayed_fac',
                     'boston60_free_fac', 'tmta_time_fac', 'tmtb_time_fac', 'kissing_dancing_images_fac',
                     'kissing_dancing_words_fac', 'action_naming_free_fac', 'g_15objects_total_np',
                     'm_wordlistwms_retention_np', 'm_wordlistwms_totalrecognition_np', 'fe_verb_fluency_np',
                     'composite_orientation', 'composite_memory_wms', 'composite_memory_fname_names', 'composite_memory_fname_occupations',
                     'composite_memory_visual', 'composite_language', 'composite_execute_functions_proc_speed',
                     'composite_visuospatial', 'composite_praxis', 'composite_attention', 'composite_global_cognition')

mlms <- function(lista_datos, datos_clinica, main_output_file = "Reporte_completo_modelos_lineales_mixtos.html", supplementary_output_file = "Resumen_reporte_modelos_lineales_mixtos.html"){
    
    resultados <- list()
    
    #datos_clinica$Plasma_AB42_40_V0_dic <- as.factor(datos_clinica$Plasma_AB42_40_V0_dic) #variable plasma como factor para ajustar el modelo
    #datos_clinica$apoe4_al <- as.factor(datos_clinica$apoe4_al) #variable apoe4_al como factor para ajustar el modelo
    
    for (i in 1:length(lista_datos)){
        
        #Definimos la variable actual dentro de la lista para cada iteracción
        variable_actual <- lista_datos[i]  
        
        #Eliminamos las filas con NA para las variables que se están empleando en la iteracción i para no perder desde el inicio tantas observaciones. En cada iteracción hay un número diferente de observaciones
        #La variable Time  y facehbi_id_fac no tienen NAs, pero es mejor añadir todas por si se cambia el dataframe de partida y tuviera
        datos_no_na <- datos_clinica[complete.cases(datos_clinica[ , c(variable_actual, "basal_age", "apoe4_al", "Plasma_AB42_40_V0_dic", "Time", "facehbi_id_fac")]), ] 
        
        # El término fixed es dinámico, cambia en cada iteracción, por lo que hay usar as.formula
        formula_modelo <- as.formula(paste(variable_actual, 
                                           "~ basal_age + apoe4_al + Plasma_AB42_40_V0_dic*Time"))
        
        # Construir el modelo
        modelo <- lme(fixed = formula_modelo,  #  fixed effects formula
                      random = ~ Time | facehbi_id_fac,            # Random intercept and slopes for subjects
                      correlation = corCompSymm(form = ~ Time | facehbi_id_fac),  # Compound symmetry structure
                      data = datos_no_na, #Dataset
                      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)) #Para optimizar los parámetros
        
        
        #Análisis de simple slopes
        simples_slopes <- emtrends(modelo, specs=c("Plasma_AB42_40_V0_dic"), var="Time")
        p_valores <- test(emtrends(modelo, specs=c("Plasma_AB42_40_V0_dic"), var="Time"))
        
        #Análisis de efectos simples
        efectos_simples <- emmeans(modelo, specs=c("Plasma_AB42_40_V0_dic", "Time"), at=list(Time=c(0,24,60)))
        efectos_simples_2 <-contrast(efectos_simples, method="pairwise", simple="Plasma_AB42_40_V0_dic")
        
        #Gráficas de predicción
        pred_plasma <- ggpredict(modelo, "Plasma_AB42_40_V0_dic") %>% plot()
        pred_time_plasma <- ggpredict(modelo, c("Time","Plasma_AB42_40_V0_dic")) %>% plot()
        pred_plasma_time <- ggpredict(modelo, c("Plasma_AB42_40_V0_dic", "Time")) %>% plot()

        # Guardar el resumen del modelo en la lista de resultados
        resultados[[i]] <- list(variable = variable_actual, 
                                datos = as.data.frame(cbind(coef(summary(modelo)),intervals(modelo,which = "fixed")$fixed[, c(1,3)])),
                                est = performance(modelo), 
                                simples_slopes = simples_slopes,
                                p_valores = p_valores,
                                efectos_simples = efectos_simples_2,
                                pred_plasma = pred_plasma, 
                                pred_time_plasma = pred_time_plasma, 
                                pred_plasma_time = pred_plasma_time,
                                #Para el html resumido
                                aic = as.data.frame(performance(modelo)$AIC), 
                                fila_plasma = as.data.frame(t(coef(summary(modelo))[6,]))
                                )
        
    }
    
    # Guardar la lista de resultados en un archivo RDS para pasar a RMarkdown
    saveRDS(resultados, file = "resultados_modelos.rds")
    
    # Renderizar el archivo RMarkdown en HTML usando los datos guardados
    rmarkdown::render("reporte_resultados_modelos.Rmd", output_file = main_output_file)
    rmarkdown::render("resumen_reporte_resultados_modelos.Rmd", output_file = supplementary_output_file)
    
    return(invisible(resultados))
    
}


# Definir el contenido del archivo .Rmd para el reporte completo

contenido_rmd_completo <- "
---
title: \"Reporte completo de Resultados de Modelos Lineales Mixtos\"
date: \"`r Sys.Date()`\"
output:
  html_document:
---

# Introduccion

Este reporte muestra los resultados de los modelos linales mixtos ajustados para las siguientes variables: *'fbb_pet_centilod', 'aHipocampal_vol', 'eTIV', 'aVentricular_vol', 'Cortical_vol', 
                     'mmse_fac', 'fname_total_fac', 'fe_skt_time_np', 'rbans_copy_fac', 'rbans_delayed_fac',
                     'boston60_free_fac', 'tmta_time_fac', 'tmtb_time_fac', 'kissing_dancing_images_fac',
                     'kissing_dancing_words_fac', 'action_naming_free_fac', 'g_15objects_total_np',
                     'm_wordlistwms_retention_np', 'm_wordlistwms_totalrecognition_np', 'fe_verb_fluency_np',
                     'composite_orientation', 'composite_memory_wms', 'composite_memory_fname_names', 'composite_memory_fname_occupations',
                     'composite_memory_visual', 'composite_language', 'composite_execute_functions_proc_speed',
                     'composite_visuospatial', 'composite_praxis', 'composite_attention', 'composite_global_cognition'*.  

Para la construcción del modelo se ha empleado la siguiente fórmula:  

```{r, eval=FALSE}
mod <- lme(
  fixed = variable_actual ~ basal_age + apoe4_al + Plasma_AB42_40_V0_dic + Time + Plasma_AB42_40_V0_dic:Time,
  random = ~ Time | facehbi_id_fac,             
  correlation = corCompSymm(form = ~ Time | facehbi_id_fac),  
  data = datos_clinica,  
  control = lmeControl(opt = \"optim\", maxIter = 100, msMaxIter = 100))
```

```{r setup, include=FALSE}

knitr::opts_chunk$set(echo = TRUE)
```

```{r results='asis', echo=FALSE}
resultados <- readRDS(\"resultados_modelos.rds\")

for (i in 1:length(resultados)) {
    
    cat(paste0(\"## Modelo para la variable: \", resultados[[i]]$variable, \"\\n\"))

    cat(knitr::kable(resultados[[i]]$datos, digits = 3, caption = \"Resumen de coeficientes del modelo\") %>%
          kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

    cat(knitr::kable(resultados[[i]]$est, digits = 3, caption = \"Estadisticas de rendimiento del modelo\") %>%
          kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

    cat(knitr::kable(resultados[[i]]$simples_slopes, digits = 3, caption = \"Analisis de pendientes simples\") %>%
          kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

    cat(knitr::kable(resultados[[i]]$p_valores, digits = 3, caption = \"p-valores de analisis de pendientes simples\") %>%
          kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

    cat(knitr::kable(resultados[[i]]$efectos_simples, digits = 3, caption = \"Analisis de efectos simples\") %>%
          kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

    print(resultados[[i]]$pred_plasma + labs(title = \"Prediccion de Plasma_AB42_40_V0_dic\"))
    
    
    print(resultados[[i]]$pred_time_plasma + labs(title = \"Prediccion de Time y Plasma_AB42_40_V0_dic\"))

    
    print(resultados[[i]]$pred_plasma_time + labs(title = \"Prediccion de Plasma_AB42_40_V0_dic y Time\"))

    cat(\"\\n---\\n\\n\") 
                               
}

```
"

# Definir el contenido del archivo .Rmd para el reporte resumido

contenido_rmd_resumen <- "
---
title: \"Resumen de los principales parámetros de Resultados de Modelos Lineales Mixtos\"
date: \"`r Sys.Date()`\"
output:
  html_document:
---

# Introduccion

Este reporte muestra el resumen de los parametros mas importantes de los modelos linales mixtos ajustados para las siguientes variables: *'fbb_pet_centilod', 'aHipocampal_vol', 'eTIV', 'aVentricular_vol', 'Cortical_vol', 
                     'mmse_fac', 'fname_total_fac', 'fe_skt_time_np', 'rbans_copy_fac', 'rbans_delayed_fac',
                     'boston60_free_fac', 'tmta_time_fac', 'tmtb_time_fac', 'kissing_dancing_images_fac',
                     'kissing_dancing_words_fac', 'action_naming_free_fac', 'g_15objects_total_np',
                     'm_wordlistwms_retention_np', 'm_wordlistwms_totalrecognition_np', 'fe_verb_fluency_np',
                     'composite_orientation', 'composite_memory_wms', 'composite_memory_fname_names', 'composite_memory_fname_occupations',
                     'composite_memory_visual', 'composite_language', 'composite_execute_functions_proc_speed',
                     'composite_visuospatial', 'composite_praxis', 'composite_attention', 'composite_global_cognition'*.  


```{r setup, include=FALSE}

knitr::opts_chunk$set(echo = TRUE)
```

```{r results='asis', echo=FALSE}

resultados <- readRDS(\"resultados_modelos.rds\")

tabla_aic <- data.frame(Variable = character(), AIC = numeric(), stringsAsFactors = FALSE)

tabla_coeficientes_plasma <- data.frame(Variable = character(), Value = numeric(), 
                                        Std.Error = numeric(), DF = numeric(), 
                                        t.value = numeric(), p.value = numeric(),
                                        stringsAsFactors = FALSE)

for (i in 1:length(resultados)) {

    variable_actual <- resultados[[i]]$variable
    
    
    colnames(resultados[[i]]$aic) <- \"AIC\" 
    tabla_aic <- rbind(tabla_aic, 
                 data.frame(Variable = variable_actual, AIC = resultados[[i]]$aic[, \"AIC\"]))

    
    fila_plasma_actual <- resultados[[i]]$fila_plasma
    colnames(fila_plasma_actual) <- c(\"Value\", \"Std.Error\", \"DF\", \"t-value\", \"p-value\")
    fila_plasma_actual$Variable <- variable_actual
    tabla_coeficientes_plasma <- rbind(tabla_coeficientes_plasma, fila_plasma_actual)
    
}

    cat(\"## Tabla de AIC por variable\n\")
    cat(knitr::kable(tabla_aic, digits = 3) %>% kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))


    cat(\"## Tabla de coeficientes para la interacción Plasma_AB42_40_V0_dic:Time\n\")
    cat(knitr::kable(tabla_coeficientes_plasma, digits = 3) %>% kable_styling(bootstrap_options = c(\"striped\", \"hover\", \"condensed\"), full_width = T))

```
"


writeLines(contenido_rmd_completo, "reporte_resultados_modelos.Rmd")
writeLines(contenido_rmd_resumen, "resumen_reporte_resultados_modelos.Rmd")

mlms(lista_variables, dfl, main_output_file = "Reporte_completo_modelos_lineales_mixtos.html", supplementary_output_file = "Resumen_reporte_modelos_lineales_mixtos.html")
