# Version simplifiée de la fonction pretreatment()
pretreatment <- function(params, serie, Int = TRUE, first = NULL, last = NULL, yrInterp = 10) {
  # params : matrice avec colonnes : top depth, bottom depth, top age, bottom age, volume
  # serie  : série brute (ex. comptages)
  
  # 1. Extraire les données
  top_depth <- params[, 1]
  bottom_depth <- params[, 2]
  top_age <- params[, 3]
  bottom_age <- params[, 4]
  volume <- params[, 5]
  
  # 2. Calcul âge moyen et épaisseur pour chaque couche
  mean_age <- (top_age + bottom_age) / 2
  thickness <- bottom_depth - top_depth
  
  # 3. Calcul des taux d'accumulation bruts
  acc <- serie / (bottom_age - top_age)   # série / durée
  
  # 4. Définir la grille d'interpolation
  if (is.null(first)) first <- min(mean_age)
  if (is.null(last))  last  <- max(mean_age)
  ages_interp <- seq(first, last, by = yrInterp)
  
  # 5. Interpolation si demandé
  if (Int) {
    acc_interp <- approx(mean_age, acc, xout = ages_interp, method = "linear", rule = 2)$y
    depth_interp <- approx(mean_age, (top_depth + bottom_depth)/2, xout = ages_interp, method = "linear", rule = 2)$y
  } else {
    acc_interp <- acc
    depth_interp <- (top_depth + bottom_depth) / 2
    ages_interp <- mean_age
  }
  
  # 6. Retourner les résultats
  return(list(
    cmI = depth_interp,  # profondeurs interpolées
    ybpI = ages_interp,  # âges interpolés
    accI = acc_interp    # taux d'accumulation interpolés
  ))
}

# --- Exemple d'utilisation ---
# Données fictives
params <- matrix(c(
  0, 1, 0, 50, 10,
  1, 2, 50, 100, 10,
  2, 3, 100, 150, 10,
  3, 4, 150, 200, 10
), ncol = 5, byrow = TRUE)

serie <- c(5, 10, 15, 20)

# Appel de la fonction
res <- pretreatment(params, serie, Int = TRUE, first = 0, last = 200, yrInterp = 10)


