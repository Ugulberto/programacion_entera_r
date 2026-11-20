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
coef <- c(11.46,13.66,14.3,12.38,14.76,11.74,-18.5,-11.2,-10.3,-5,-5.1)
A <- rbind(c(1,rep(0,10)),
           c(0,1,rep(0,9)),
           c(0,0,1,rep(0,8)),
           c(0,0,0,1,rep(0,7)),
           c(0,0,0,0,1,rep(0,6)),
           c(0,0,0,0,0,1,rep(0,5)),
           c((rep(1,6)),rep(0,5)),
           c(0,-1,rep(0,4),29,rep(0,4)),
           c(0,0,-1,0,0,0,0,29,0,0,0),
           c(0,0,0,-1,0,0,0,0,29,0,0),
           c(0,0,0,0,-1,0,0,0,0,29,0),
           c(0,0,0,0,0,-1,0,0,0,0,29))
b <- c(35,30,40,50,35,25,180,0,0,0,0,0)
dir <- c("<=","<=","<=","<=","<=","<=","=",">=",">=",">=",">=",">=")

resolver_programacion_entera(type, coef, A, dir, b)

cat("\n\n--- MEJOR SOLUCION ---\n")
imprimir_solucion(maxSol, maxVars)