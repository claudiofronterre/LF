# Load required packages and functions ----------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("GGally", "ggcorrplot", "dplyr", "DiagrammeR", "corrplot", "tidyverse", "mgcv") # package names
pacman::p_load(pkgs, character.only = T)

# Load the data -----------------------------
onco <- read.csv("data/Ethiopia_ONC_CommunityData.csv")
covariates <- readRDS("data/covariates.rds")
covariates[618, 41] <- 10.18987 # this values has been obtained averaging the neighboring pixels

# Create new urbanisation index  and population year based on the survey year ------------------
urban <- covariates[, c("Urbanization2000", "Urbanization2010", "Urbanization2015")]
pop <- covariates[, c("PopDens2000", "PopDens10adj", "PopDens15adj")]
id <- as.factor(onco$surveyyear)
levels(id) <- c("2000", "2010", "2010", "2015")
id <- as.numeric(id)
covariates <- covariates[ , -grep(pattern = "Urb", x = names(covariates))]
covariates <- covariates[ , -grep(pattern = "Pop", x = names(covariates))]
covariates <- covariates[ , -which(names(covariates) == "GLC2006")]
covariates$Urbanisation <- apply(cbind(urban, id), 1, function (x) x[x["id"]])
covariates$PopDens <- apply(cbind(pop, id), 1, function (x) x[x["id"]])

# Exploratory analysis and multicollinearty check ------------------

# Grouping the variables by category (climate, demo, hydro, topo, veg)
temp <- c(2:11, 19:21)
precip <- c(1, 12:18)
climate <- 1:21
demo <- c(22, 37, 38)
hydro <- 23:27
topo <- 28:30
veg <- 31:36
groups <- list(Climate = climate, Demographic = demo, Hydrology = hydro, Topographic = topo, Vegetation = veg)
groups <- list(Temperature = temp, Demographic = demo, 
               Precipitation = precip, Hydrology = hydro, 
               Topographic = topo, Vegetation = veg)


# Convert the non numeric variables to factors
covariates$Drainage_Basins <- as.factor(covariates$Drainage_Basins)
covariates$GLC2006_reclass <- as.factor(covariates$GLC2006_reclass)
covariates$Urbanisation <- as.factor(covariates$Urbanisation)
saveRDS(covariates, "data/covariates_final.rds")

# Create ids for the columns that contains numeric and categorical variables
id_num <- as.numeric(which(sapply(covariates, class) == "numeric"))
id_fac <- as.numeric(which(sapply(covariates, class) != "numeric"))


#Transform some of the covariates
#covariates$EucDist_km_Rivers <- log10(covariates$EucDist_km_Rivers + 1)

# Box plot and histograms(asses outliers, reasonable values?)
fill <- "#4271AE"
line <- "#1F3552"
for (i in 1:length(groups)) {
  df <- tidyr::gather(data.frame(scale(covariates[ , intersect(groups[[i]], id_num)]))) # scale them
  p <- ggplot(data = df, mapping = aes(x = key, y = value)) +
    geom_boxplot(fill = fill, colour = line, alpha = 0.5,
                 outlier.colour = "#1F3552", outlier.shape = 20) +
    #facet_wrap(~ key) + 
    ggtitle(paste(names(groups[i]), "variables")) +
    theme_bw()
  print(p)
}

# Correlation plot
ggcorrplot(corr = cor(covariates[,id_num]), 
           type = "lower", 
           ggtheme = ggplot2::theme_minimal,
           hc.order = T, 
           show.diag = F,
           p.mat = cor_pmat(covariates[,id_num]),
           outline.col = "white", 
           lab = F, 
           legend.title = "Correlation",
           tl.cex = 11, tl.srt = 55)
# Diagramm of correlation 
#library("corrr")
#network_plot(correlate(covariates), min_cor = 0.5, legend = F)


# Selection with GAM --------------------------
id <- onco$surveytype == "REMO"
onco <- onco[id, ]
covariates <- covariates[id, ]
onco$prevalence <- onco$NP_pos/onco$examined
df_model <- cbind(onco, covariates)

## Non spatial models -------------------------

# Pairs plot with output
for (i in 1:length(groups)) {
  df <- data.frame(prevalence = onco$prevalence, covariates[ , groups[[i]]])
  p <- ggpairs(data = df, 
          title = paste(names(groups[i]), "variables"), 
          lower = list(continuous = wrap("smooth_loess", alpha = 0.1, color = "#001782"))) + 
    theme_bw()
  print(p)
}


# Plot univariate models with every covariates
for (i in 1:length(groups)) {
  df <- data.frame(prevalence = onco$prevalence, examined = onco$examined, 
                   covariates[ , intersect(groups[[i]], id_num)])
  df <- tidyr::gather(df, key = covariate, value = value, -c(prevalence, examined))
  p <- ggplot(data = df, mapping = aes(x = value, y = prevalence)) +
    geom_point(color = fill, alpha = 0.1) +
    geom_smooth(mapping = aes(weight = examined), color = line, fill = fill, alpha = 0.5, size = 0.8, 
                method = "glm", se = T, formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit"))) +
    facet_wrap(~ covariate, scales = "free_x") + 
    ggtitle(paste(names(groups[i]), "variables")) +
    theme_bw()
  print(p)
}

# Model selection grouping the variables
models <- lapply(groups, function (x) {
  if(any(sapply(covariates[,x], class) == "factor")) {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, names(covariates)[intersect(x, id_num)], sep = "(", collapse = ")+"), ")")
    #cov_char_num <- paste(names(covariates)[intersect(x, id_num)], collapse= "+") # no spline
    cov_char_fac <- paste(names(covariates)[intersect(x, id_fac)], collapse= "+")
    f <- as.formula(paste("prevalence", paste(cov_char_num, cov_char_fac, sep = "+"), sep = "~"))
    #f <- as.formula(paste("prevalence", paste(names(covariates)[x], collapse= "+"), sep = "~")) #no spline
  } else {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, names(covariates)[intersect(x, id_num)], sep = "(", collapse = ")+"), ")")
    f <- as.formula(paste("prevalence", cov_char_num, sep = "~"))
    #f <- as.formula(paste("prevalence", paste(names(covariates)[x], collapse= "+"), sep = "~")) #no spline
  }
  gam(formula = f, weights = examined, data = df_model, family = quasibinomial(), 
      method = "REML", select = T, gamma = 1.5)
})

lapply(models, summary)
lapply(models, anova)
Aic <- sapply(models, function(x) x$aic)
dev_explained <- sapply(models, function(x) summary(x)$dev.expl)
pvalue <- sapply(models, function(x) round(as.numeric(summary(x)$p.pv[2]), 4))


# Add interaction to the formula
#int_basis <- "ti"
#f_int <- as.name(paste(combn(length(names(covariates)[veg]), m = 2, FUN = function(x) {
#  paste0(int_basis, "(", paste(names(covariates)[veg][x], collapse = ","), ")")
#}), collapse = "+"))

# par(mfrow = c(2, ncol(covariates[, intersect(g, id_num)])))
# plot(fit2, n = 1000) # n increased to make plot smooth
# plot(fit2, residuals = TRUE, n = 1000)

# No grouping
all_cov <- as.list(names(covariates)) # put all the names of the covariates in a list

models <- lapply(all_cov, function (x) {
  if(class(covariates[,x]) == "numeric") {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, x, sep = "("), ")")
    f <- as.formula(paste("prevalence", cov_char_num, sep = "~"))
    #f <- as.formula(paste("prevalence", x, sep = "~")) # no spline
  } else {
    f <- as.formula(paste("prevalence", x, sep = "~"))
  }
  gam(formula = f, weights = examined, data = df_model, family = quasibinomial(), method="REML", select = T)
})

lapply(models, function(x) try(plot(x)))
sapply(models, gam.check)
Aic <- sapply(models, function(x) x$aic)
dev_explained <- sapply(models, function(x) summary(x)$dev.expl)
pvalue <- sapply(models, function(x) round(as.numeric(summary(x)$p.pv[2]), 4))
r2 <- sapply(models, function(x) summary(x)$r.sq)

group_names <- rep(0, ncol(covariates))
for (i in 1:length(groups)) group_names[groups[[i]]] <- names(groups[i])

results <- data.frame(group = group_names, covariates = names(covariates), aic = Aic, dev_expl = dev_explained, pvalue, r2)

res_plot <- results[order(results$dev_expl, decreasing = F), ]
res_plot$covariates <- factor(res_plot$covariates, levels = res_plot$covariates)


# set ggthemr theme
#ggthemr::ggthemr("flat", type = "outer", layout = "scientific", spacing = 2) 
varImp <- ggplot(res_plot, aes(x = covariates, y = dev_expl)) + 
            geom_segment(aes(x = covariates, 
                             xend = covariates, 
                             y = 0, 
                             yend = dev_expl), 
                             color = "black",
                             linetype = "solid", 
                             size = 0.3) +   # Draw dashed lines
            geom_point(aes(fill = group), color = "black", shape = 21, size = 3) +   # Draw points
            labs(title="Variable Importance", 
                 subtitle="Ordered by Deviance Explained") +
            coord_flip() + 
            theme_classic(base_size = 13)
varImp + ggsci::scale_fill_locuszoom(name = "")
# remove all ggthemr effects 
#ggthemr::ggthemr_reset()

selected_variables <- as.character(results$covariates[order(results$dev_expl, decreasing = T)])
saveRDS(selected_variables, file = "output/selected_variables.rds")

# Spatial models

# No grouping
models <- lapply(all_cov, function (x) {
  if(class(covariates[,x]) == "numeric") {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, x, sep = "("), ")")
    f <- as.formula(paste("prevalence", cov_char_num, sep = "~"))
    f <- update(f, ~ . + s(longitude, latitude))
  } else {
    f <- as.formula(paste("prevalence", x, sep = "~"))
    f <- update(f, ~ . + s(longitude, latitude))
  }
  gam(formula = f, weights = examined, data = df_model, family = quasibinomial(), method="REML", select = T)
})

lapply(models, function(x) try(plot(x)))
lapply(models, summary)
lapply(models, anova)
par(mfrow = c(2, 2))
sapply(models, gam.check)
Aic <- sapply(models, function(x) x$aic)
dev_explained <- sapply(models, function(x) summary(x)$dev.expl)
pvalue <- sapply(models, function(x) round(as.numeric(summary(x)$p.pv[2]), 4))
r2 <- sapply(models, function(x) summary(x)$r.sq)

results <- data.frame(group = group_names, covariates = names(covariates), aic = Aic, dev_expl = dev_explained, r2)

res_plot <- results[order(results$dev_expl, decreasing = F), ]
res_plot$covariates <- factor(res_plot$covariates, levels = res_plot$covariates)


# set ggthemr theme
ggthemr::ggthemr("flat", type = "outer", layout = "scientific", spacing = 2) 
varImp_sp <- ggplot(res_plot, aes(x = covariates, y = dev_expl)) + 
  geom_segment(aes(x = covariates, 
                   xend = covariates, 
                   y = min(dev_expl), 
                   yend = dev_expl), 
                   col = "black",
                   linetype = "solid", 
                   size = 0.3) +   # Draw dashed lines
  geom_point(aes(color = group), size = 3) +   # Draw points
  labs(title="Variable Importance", 
       subtitle="Ordered by Deviance Explained") +
  coord_flip() + 
  theme_classic()
varImp_sp
# remove all ggthemr effects 
ggthemr::ggthemr_reset()




# Model selection grouping the variables
models <- lapply(groups, function (x) {
  if(any(sapply(covariates[,x], class) == "factor")) {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, names(covariates)[intersect(x, id_num)], sep = "(", collapse = ")+"), ")")
    #cov_char_num <- paste(names(covariates)[intersect(x, id_num)], collapse= "+") # no spline
    cov_char_fac <- paste(names(covariates)[intersect(x, id_fac)], collapse= "+")
    f <- as.formula(paste("prevalence", paste(cov_char_num, cov_char_fac, sep = "+"), sep = "~"))
    #f <- as.formula(paste("prevalence", paste(names(covariates)[x], collapse= "+"), sep = "~")) #no spline
    f <- update(f, ~ . + s(longitude, latitude))
  } else {
    basis <- "s"
    cov_char_num <- paste0(paste(basis, names(covariates)[intersect(x, id_num)], sep = "(", collapse = ")+"), ")")
    f <- as.formula(paste("prevalence", cov_char_num, sep = "~"))
    #f <- as.formula(paste("prevalence", paste(names(covariates)[x], collapse= "+"), sep = "~")) #no spline
    f <- update(f, ~ . + s(longitude, latitude))
  }
  gam(formula = f, data = df_model, weights = examined, family = quasibinomial(),
      method = "REML", select = F)
})

#x11()
#plot(models, pages = 1, scale = 0)
par(mfrow = c(2, 2))
lapply(models, gam.check)
sapply(models, function(x) try(plot(x, scale = 0, page = 1)))
lapply(models, summary)
lapply(models, anova)
Aic <- sapply(models, function(x) x$aic)
dev_explained <- sapply(models, function(x) summary(x)$dev.expl)
pvalue <- sapply(models, function(x) round(as.numeric(summary(x)$s.pv), 4))
pvalue <- sapply(pvalue, function(x) x[-length(x)])
pstar <- sapply(pvalue, function(x) which(x < .15))
groups

id_selected <- purrr::map2(groups, pstar, ~ .x[.y]) %>%
  purrr::reduce(.f = c)
names(covariates)[id_selected]

covariates_final <- covariates[, id_selected]
cor(covariates_final)

ggcorrplot(corr = cor(covariates_final), 
           type = "lower", 
           ggtheme = ggplot2::theme_minimal,
           hc.order = T, 
           show.diag = F,
           p.mat = cor_pmat(covariates_final),
           outline.col = "white", 
           lab = T, 
           legend.title = "Correlation",
           tl.cex = 11, tl.srt = 55)
cor_mat <- cor(covariates_final)
diag(cor_mat) <- 0
idx <- which(abs(cor_mat) > .7, arr.ind = T)
col <- unique(idx[, 2])
selected_covariates <- covariates_final[, -col]

ggcorrplot(corr = cor(selected_covariates), 
           type = "lower", 
           ggtheme = ggplot2::theme_minimal,
           hc.order = T, 
           show.diag = F,
           p.mat = cor_pmat(selected_covariates),
           outline.col = "white", 
           lab = T, 
           legend.title = "Correlation",
           tl.cex = 11, tl.srt = 55)

saveRDS(names(selected_covariates), file = "output/selected_variables.rds")
