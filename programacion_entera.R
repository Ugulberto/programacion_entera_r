# Creado por: David Blanco
# Método de programación lineal entera branch and bound
# Se recomienda ejecutar en modo no echo (Ctrl + A y luego Ctrl + Shift + S)

install.packages("lpSolve")
library(lpSolve)

max <- -Inf
maxVars <- NA
maxSol <- NA

# Helper para repetir la sangría según la profundidad del árbol
indentar <- function(nivel) {
  if (nivel > 0) cat(rep("|   ", nivel), sep = "")
}

imprimir_problema <- function(tipo_optimizacion, coeficientes, A, desigualdades, b, nivel = 0) {
  indentar(nivel); cat("Problema actual (Nivel de profundidad de rama ", nivel, "):\n", sep = "")
  indentar(nivel); cat(tipo_optimizacion, "imizar z = ", sep = "")
  
  for (i in seq_along(coeficientes)) {
    if (i > 1) cat(" + ")
    cat(coeficientes[i], "x_", i, sep = "")
  }
  
  cat("\n")
  indentar(nivel); cat("Sujeto a:\n")
  for (i in 1:nrow(A)) {
    indentar(nivel)
    primer_termino <- TRUE
    for (j in 1:ncol(A)) {
      if (A[i, j] != 0) {
        if (!primer_termino) cat(" + ")
        cat(A[i, j], "x_", j, sep = "")
        primer_termino <- FALSE
      }
    }
    cat(" ", desigualdades[i], " ", b[i], "\n", sep = "")
  }
}

imprimir_solucion <- function(solucion, variables, nivel = 0) {
  indentar(nivel); cat("Solución: z = ", solucion, " con las variables: ", sep = "")
  for (i in seq_along(variables)) {
    if (i > 1) cat(", ")
    cat("x_", i, " = ", round(variables[i], 4), sep = "")
  }
  cat("\n")
}

resolver_programacion_entera <- function(tipo_optimizacion, coeficientes, A, desigualdades, b, nivel = 0) {
  cat("\n")
  imprimir_problema(tipo_optimizacion, coeficientes, A, desigualdades, b, nivel)
  solucion <- lp(tipo_optimizacion, coeficientes, A, desigualdades, b)
  
  if (solucion$status == 2) {
    indentar(nivel); cat("[X] Problema infactible (Fin de rama).\n")
  } else {
    # Problema factible.
    if (all(abs(solucion$solution - round(solucion$solution)) < 1e-6)) {
      # Solución válida
      imprimir_solucion(solucion$objval, solucion$solution, nivel)
      if (is.nan(max) || solucion$objval > max) {
        indentar(nivel); cat("[*] Nuevo mejor máximo entero: z = ", solucion$objval, "\n", sep = "")
        max <<- solucion$objval
        maxVars <<- solucion$solution
        maxSol <<- solucion$objval
      }
      indentar(nivel); cat("-------------------------------------\n")
    } else {
      # Solución con números decimales en las variables
      imprimir_solucion(solucion$objval, solucion$solution, nivel)
      if (solucion$objval < max) {
        indentar(nivel); cat("[-] Poda por cota: tenemos un máximo mejor (", max, "), paramos rama.\n", sep = "")
        indentar(nivel); cat("-------------------------------------\n")
        return()
      }
      
      minDist <- Inf
      minVarInd <- -1
      
      for (i in seq_along(solucion$solution)) {
        num <- solucion$solution[i]
        if (abs(num - round(num)) > 1e-6) {
          if (abs(floor(num) + 0.5 - num) < minDist) {
            minDist <- abs(floor(num) + 0.5 - num)
            minVarInd <- i
          }
        }
      }
      
      num <- solucion$solution[minVarInd]
      
      # Rama izquierda (<=)
      indentar(nivel); cat("--- Rama 1: Añadimos restricción x_", minVarInd, " <= ", floor(num), "\n", sep = "")
      resolver_programacion_entera(tipo_optimizacion, 
                                   coeficientes,
                                   rbind(A, c(rep(0, minVarInd - 1), 1, rep(0, length(coeficientes) - minVarInd))),
                                   c(desigualdades, "<="),
                                   c(b, floor(num)),
                                   nivel + 1)
      
      # Rama derecha (>=)
      indentar(nivel); cat("+-- Rama 2: Añadimos restricción x_", minVarInd, " >= ", ceiling(num), "\n", sep = "")
      resolver_programacion_entera(tipo_optimizacion, 
                                   coeficientes,
                                   rbind(A, c(rep(0, minVarInd - 1), 1, rep(0, length(coeficientes) - minVarInd))),
                                   c(desigualdades, ">="),
                                   c(b, ceiling(num)),
                                   nivel + 1)
    }
  }
}

type <- "max"
coef <- c(11.46, 13.66, 14.3, 12.38, 14.76, 11.74, # Beneficio por vender una caja en el pueblo i
          -9.46, -14.77, -4.06, -2.81, -6.29,      # Coste de realizar la ruta r una vez
          0, 0, 0, 0, 0, 0, 0)                     # q_ir (qA1, qA2, qB2, qC3, qD4, qD5, qE5), no aportan beneficio directo

A <- rbind(
  c(1,0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el origen O (<= 35)
  c(0,1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el pueblo A (<= 30)
  c(0,0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el pueblo B (<= 40)
  c(0,0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el pueblo C (<= 50)
  c(0,0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el pueblo D (<= 35)
  c(0,0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,0,0), # Límite de cajas en el pueblo E (<= 25)
  
  c(1,1,1,1,1,1,  0,0,0,0,0,  0,0,0,0,0,0,0), # Se deben vender las 180 cajas en total
  
  c(0,0,0,0,0,0, -29,0,0,0,0,  1,0,0,0,0,0,0), # Capacidad de la ruta 1 (qA1 - 29*y1 <= 0)
  c(0,0,0,0,0,0,  0,-29,0,0,0, 0,1,1,0,0,0,0), # Capacidad de la ruta 2 (qA2 + qB2 - 29*y2 <= 0)
  c(0,0,0,0,0,0,  0,0,-29,0,0, 0,0,0,1,0,0,0), # Capacidad de la ruta 3 (qC3 - 29*y3 <= 0)
  c(0,0,0,0,0,0,  0,0,0,-29,0, 0,0,0,0,1,0,0), # Capacidad de la ruta 4 (qD4 - 29*y4 <= 0)
  c(0,0,0,0,0,0,  0,0,0,0,-29, 0,0,0,0,0,1,1), # Capacidad de la ruta 5 (qD5 + qE5 - 29*y5 <= 0)
  
  c(0,-1,0,0,0,0, 0,0,0,0,0,  1,1,0,0,0,0,0), # Cajas vendidas en A mediante rutas 1 y 2 (qA1 + qA2 - xA = 0)
  c(0,0,-1,0,0,0, 0,0,0,0,0,  0,0,1,0,0,0,0), # Cajas vendidas en B mediante ruta 2 (qB2 - xB = 0)
  c(0,0,0,-1,0,0, 0,0,0,0,0,  0,0,0,1,0,0,0), # Cajas vendidas en C mediante ruta 3 (qC3 - xC = 0)
  c(0,0,0,0,-1,0, 0,0,0,0,0,  0,0,0,0,1,1,0), # Cajas vendidas en D mediante rutas 4 y 5 (qD4 + qD5 - xD = 0)
  c(0,0,0,0,0,-1, 0,0,0,0,0,  0,0,0,0,0,0,1)  # Cajas vendidas en E mediante ruta 5 (qE5 - xE = 0)
)

b <- c(35, 30, 40, 50, 35, 25, # Límites de almacén
       180,                    # Total de cajas
       0, 0, 0, 0, 0,          # Capacidad de las 5 rutas
       0, 0, 0, 0, 0)          # Enlace q_ir - x_i = 0

dir <- c("<=", "<=", "<=", "<=", "<=", "<=",
         "=",
         "<=", "<=", "<=", "<=", "<=",
         "=", "=", "=", "=", "=")

resolver_programacion_entera(type, coef, A, dir, b)

cat("\n\n--- MEJOR SOLUCION ---\n")
imprimir_solucion(maxSol, maxVars)
