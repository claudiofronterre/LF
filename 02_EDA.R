# Load required packages and functions -------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("ggplot2") # package names
pacman::p_load(pkgs, character.only = T)

# Load data ---------------------------------------------
lf_both <- readRDS("data/lf_both.rds")

# Exploratory analysis ----------------------------------
p_ict <- lf_both$Prevalence_ict
p_mf <- lf_both$Prevalence_mf
elog_ict <- log((lf_both$Positive_ict + 0.5)/(lf_both$Examined - lf_both$Positive_ict + 0.5))
elog_mf <- log((lf_both$Positive_mf + 0.5)/(lf_both$Examined - lf_both$Positive_mf + 0.5))

# Tables with summary statistics
quant <- rbind(quantile(elog_ict, probs = seq(0, 1, 0.25)), quantile(elog_mf, probs = seq(0, 1, 0.25)))
tabl <- data.frame(Prevalence = c("ICT", "MF"), Mean = c(mean(elog_ict), mean(elog_mf)), Sd = c(sd(elog_ict), sd(elog_mf)), quant)
names(tabl)[4:8] <- c("Min", "1Q", "Median", "3Q", "Max")
knitr::kable(tabl, digits = 3, caption = "\\label{tab:summary}Summary statistics for MF and NP prevalence on the empirical logit scale")
quant <- rbind(quantile(p_ict, probs = seq(0, 1, 0.25)), quantile(p_mf, probs = seq(0, 1, 0.25)))
tabl <- data.frame(Prevalence = c("NP", "MF"), Mean = c(mean(p_ict), mean(p_mf)), Sd = c(sd(p_ict), sd(p_mf)), quant)
names(tabl)[4:8] <- c("Min", "1Q", "Median", "3Q", "Max")
knitr::kable(tabl, digits = 3, caption = "\\label{tab:summary2}Summary statistics for MF and NP prevalence")

# Visualise the relationship
lf_both$elog_ict <- elog_ict
lf_both$elog_mf <- elog_mf

ggplot(data = lf_both, mapping = aes(x = Prevalence_ict*100, y = Prevalence_mf*100, fill = Country)) +
  geom_abline(slope = 1, intercept = 0, col = "grey", linetype = 2) + 
  geom_point(shape = 21, col = "black", size = 2) +
  xlim(0, 100) +
  ylim(0, 100) +
  xlab("ICT prevalence (%)") +
  ylab("MF prevalece (%)") + 
  scale_fill_brewer(type = "qualitative", palette = "Set3") + 
  theme_bw() +
  theme(legend.position = "top")

r1 <- ggplot(data = lf_both, mapping = aes(x = elog_ict, y = elog_mf, fill = Country)) +
  geom_smooth(aes(x = elog_ict, y = elog_mf), inherit.aes = F, method = "lm", col = "black") +
  geom_abline(slope = 1, intercept = 0, col = "grey", linetype = 2) + 
  geom_point(shape = 21, col = "black", size = 2) +
  coord_fixed(xlim = c(-7, 1.5), ylim = c(-7, 1.5)) +
  xlab("ICT empirical logit") +
  ylab("MF empirical logit") + 
  scale_fill_brewer(type = "qualitative", palette = "Set3") + 
  theme_bw() +
  theme(legend.position = "top")
r1

ggsave(filename = "figs/empirical_relationship_ict_mf.pdf", device = "pdf", width = 7, height = 7)

r2 <- ggplot(data = lf_both, mapping = aes(x = elog_ict, y = elog_mf, fill = Country)) +
  geom_smooth(aes(x = elog_ict, y = elog_mf, fill = Country), col = "black", inherit.aes = F, method = "lm") +
  geom_abline(slope = 1, intercept = 0, col = "black", linetype = 2) + 
  geom_point(shape = 21, col = "black", size = 2, alpha = 1) +
  coord_fixed(xlim = c(-7, 1.5), ylim = c(-7, 1.5)) +
  xlab("ICT empirical logit") +
  ylab("MF empirical logit") + 
  facet_wrap(~ Country) +
  scale_fill_brewer(type = "qualitative", palette = "Set3") + 
  scale_color_brewer(type = "qualitative", palette = "Set3") + 
  theme_bw() +
  guides(fill = "none", colour = "none") + 
  theme(legend.position = "top")
r2

ggsave(filename = "figs/empirical_relationship_ict_mf_byCountry.pdf", device = "pdf", width = 7, height = 7)

ggpubr::ggarrange(r1, r2 + ggpubr::rremove("xlab") + ggpubr::rremove("ylab"), common.legend = F)

ggsave(filename = "figs/empirical_relationship_ict_mf_combined.pdf", device = "pdf", width = 15, height = 8)



# Binomial
lf_both$group <- 1:nrow(lf_both)

ict_f <- bf(Positive_ict | trials(Examined) ~ (1|p|group))
mf_f <- bf(Positive_mf | trials(Examined) ~ (1|p|group))
f <- mvbf(ict_f, mf_f)

fit <- brm(formula = f, family = binomial, data = lf_both, chains = 8, cores = 8, iter = 4000)
summary(fit, R2 = T)
pp_check(fit, resp = "Positivemf")


#Intercept only model
ict_f <- bf(elog_ict ~ 1)
mf_f <- bf(elog_mf ~ 1)
f <- mvbf(ict_f, mf_f)

fit <- brm(formula = f, data = lf_both, chains = 8, cores = 8)
summary(fit)
pp_check(fit, resp = "elogict", nsamples = 50)
pp_check(fit, resp = "elogmf", nsamples = 50)
gridExtra::grid.arrange(pp1, pp2, ncol = 2)
bayes_R2(fit)

ict_f <- bf(elog_ict ~ (1|p|Country))
mf_f <- bf(elog_mf ~ (1|p|Country))
f <- mvbf(ict_f, mf_f)

fit <- brm(formula = f, data = lf_both, chains = 8, cores = 8)
summary(fit)

# Check convergence
my_sso <- launch_shinystan(fit)
deploy_shinystan(sso = my_sso, appName = "onco")



