library(ggplot2)
library(readxl)
library(ggpubr)
library(CHNOSZ)
library(dplyr)
library(ggpmisc)
library(scales)

theme_set(theme_classic())

setwd("~/Desktop/CLB23/Data")

#importing data
CLB23_sulfate <- read_excel("CLB23_data.xlsx", sheet = 8)
CLB23_sulfide <- read_excel("CLB23_data.xlsx", sheet = 9)
CLB23_methane <- read_excel("CLB23_data.xlsx", sheet = 7)
CLB23_dic <- read_excel("CLB23_data.xlsx", sheet = 6)
CLB23_porosity <- read_excel("CLB23_data.xlsx", sheet = 10)
CLB23_hydrogen <- read_excel("CLB23_data.xlsx", sheet = 3)
CLB23_cellcounts <- read_excel("CLB23_data.xlsx", sheet = 12, range = "A11:L82")
CLB23_sulfate_dic <- read_excel("CLB23_data.xlsx", sheet = 16)
CLB23_deltaG_methanogenesis <- read_excel("CLB23_data.xlsx", sheet = 4)
CLB23_deltaG_sulfatereduction <- read_excel("CLB23_data.xlsx", sheet = 5)
WOR19 <- read_excel("CLB23_data.xlsx", sheet = 14)
WOR13 <- read_excel("CLB23_data.xlsx", sheet = 17)
CLB23_d13CH4 <- read_excel("CLB23_data.xlsx", sheet = 18)
H2_leaktest <- read_excel("CLB23_data.xlsx", sheet = 19)

##############################
#        geochemistry        #
##############################

#t-test for hydrogen concentration difference between WOR and CLB
CLB23_hydrogen$Site <- "CLB"
CLB23_hydrogen_comparison <- CLB23_hydrogen %>% 
  select("Depth", "H2_aq_nM", "Site")
CLB23_hydrogen_comparison <- subset.data.frame(CLB23_hydrogen_comparison, CLB23_hydrogen_comparison$Depth < 15)
WOR19$Site <- "WOR"
WOR19_hydrogen_comparison <- WOR19 %>% 
  select("Depth", "H2_nM", "Site")
WOR19_hydrogen_comparison <- subset.data.frame(WOR19_hydrogen_comparison, WOR19_hydrogen_comparison$Depth < 15)
WOR19_hydrogen_comparison <- rename(WOR19_hydrogen_comparison, H2_aq_nM = H2_nM)
h2_comparison <- rbind(CLB23_hydrogen_comparison, WOR19_hydrogen_comparison)

t.test(H2_aq_nM ~ Site, paired = FALSE, data = h2_comparison)

#t-test for hydrogen concentration difference between CLB 0-10 cm and CLB 10+ cm
CLB23_hydrogen$Site <- "CLB"
CLB23_hydrogen_comparison <- CLB23_hydrogen %>% 
  select("Depth", "H2_aq_nM", "Site")
CLB23_hydrogen_comparison_0_10 <- subset.data.frame(CLB23_hydrogen_comparison, CLB23_hydrogen_comparison$Depth < 10)
CLB23_hydrogen_comparison_0_10$section <- "shallow"
CLB23_hydrogen_comparison_10_52 <- subset.data.frame(CLB23_hydrogen_comparison, CLB23_hydrogen_comparison$Depth > 10)
CLB23_hydrogen_comparison_10_52$section <- "deep"

h2_comparison_CLB_sections <- rbind(CLB23_hydrogen_comparison_0_10, CLB23_hydrogen_comparison_10_52)

t.test(H2_aq_nM ~ section, paired = FALSE, data = h2_comparison_CLB_sections)

#DIC vs sulfide or sulfate to see if SR is linked to methane or DIC
CLB23_sulfate_dic <- subset.data.frame(CLB23_sulfate_dic, CLB23_sulfate_dic$Depth > 12)
CLB23_sulfate_dic <- subset.data.frame(CLB23_sulfate_dic, CLB23_sulfate_dic$Depth < 40)
CLB23_sulfate_dic <- subset.data.frame(CLB23_sulfate_dic, CLB23_sulfate_dic$Core == "CLB2023 C3")

CLB23_SRsource_plot <- ggplot(CLB23_sulfate_dic, aes(x = Sulfate_used_mM, y = DIC_used_mM)) +
  stat_poly_line(se = FALSE, linetype = "dashed", na.rm = TRUE, formula = y ~ x + 0) +
  geom_abline(intercept = 0, slope = 2, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotdash") +
  stat_poly_eq(use_label(c("eq","R2","adj.R2")), na.rm = TRUE, formula = y ~ x + 0) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core), Sulfate_used_mM, DIC_used_mM), size = 2.5) +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "∆ DIC (mM)", x = "∆ Sulfate (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8)) +
  expand_limits(y = 0, x = 20) +
  theme(legend.position = "top")

WOR13_SRsource_plot <- ggplot(WOR13, aes(x = Sulfate_used, y = DIC_used)) +
  stat_poly_line(se = FALSE, linetype = "dashed", na.rm = TRUE, formula = y ~ x + 0) +
  geom_abline(intercept = 0, slope = 2, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotdash") +
  stat_poly_eq(vjust = 12, hjust = -.3, use_label(c("eq","R2","adj.R2")), na.rm = TRUE, formula = y ~ x + 0) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core), Sulfate_used, DIC_used), size = 2.5) +
  scale_shape_manual(values = c(15,20)) +
  scale_color_manual(values = c("#FFC425", "#F8766D")) +
  labs(y = "∆ DIC (mM)", x = "∆ Sulfate (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8)) +
  expand_limits(y = 0, x = 20) +
  theme(legend.position = "top")

CLB23WOR13_SRsource_plot <- ggplot() +
  stat_poly_line(data = CLB23_sulfate_dic, se = FALSE, linetype = "dashed", na.rm = TRUE, mapping = aes(x = Sulfate_used_mM, y = DIC_used_mM), formula = y ~ x + 0) +
  geom_abline(intercept = 0, slope = 2, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotdash") +
  stat_poly_eq(data = CLB23_sulfate_dic, vjust = 0.1, mapping = aes(x = Sulfate_used_mM, y = DIC_used_mM, label = after_stat(paste(eq.label, "*\",\"*~CLB"))), na.rm = TRUE, formula = y ~ x + 0) +
  geom_point(data = CLB23_sulfate_dic, aes(shape = as.character(Core), color = as.character(Core), x = Sulfate_used_mM, y = DIC_used_mM), size = 2.5) +
  stat_poly_line(data = WOR13, color = "darkgreen", se = FALSE, linetype = "dashed", na.rm = TRUE, mapping = aes(x = Sulfate_used, y = DIC_used), formula = y ~ x + 0) +
  stat_poly_eq(data = WOR13, mapping = aes(x = Sulfate_used, y = DIC_used, label = after_stat(paste(eq.label, "*\",\"*~WOR"))), na.rm = TRUE, formula = y ~ x + 0) +
  geom_point(data = WOR13, aes(shape = as.character(Core), color = as.character(Core), x = Sulfate_used, y = DIC_used), size = 2.5) +
  scale_shape_manual(values = c(16,17,18,15,20)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF", "#FFC425", "#F8766D")) +
  labs(y = "∆ DIC (mM)", x = "∆ Sulfate (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8)) +
  expand_limits(y = 0, x = 20)

setwd("~/Desktop/CLB23/Figures")

tiff("CLB23_SRsource_test.tiff", units="in", width=4, height=4, res=300)
CLB23_SRsource_plot
dev.off()

tiff("WOR13_SRsource_test.tiff", units="in", width=4, height=4, res=300)
WOR13_SRsource_plot
dev.off()

tiff("CLB23WOR13_SRsource_test.tiff", units="in", width=6, height=4, res=300)
CLB23WOR13_SRsource_plot
dev.off()

setwd("~/Desktop/CLB23/Data")
  
#creating combined geochemistry plots
CLB23_sulfate_plot <- ggplot(data = CLB23_sulfate, aes(x = Sulfate_mM, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(y = "Depth (cm)", x = "Sulfate (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8)) 

CLB23_sulfide_plot <- ggplot(data = CLB23_sulfide, aes(x = Sulfide_µM, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "Depth (cm)", x = "Sulfide (µM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_methane_plot <- ggplot(data = CLB23_methane, aes(x = Methane_mM, y = Depth)) +
  geom_errorbar(aes(xmin = Methane_mM - Methane_stdev_mM, xmax = Methane_mM + Methane_stdev_mM), color = "#77797a") +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(y = "Depth (cm)", x = "Methane (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_d13CH4_plot <- ggplot(data = CLB23_d13CH4, aes(x = corr_iCH4, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "Depth (cm)", x = bquote(~δ^13~CH[4]), shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_d13CO2_plot <- ggplot(data = CLB23_d13CH4, aes(x = corr_iCO2, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "Depth (cm)", x = bquote(~δ^13~CO[2]), shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_dic_plot <- ggplot(data = subset.data.frame(CLB23_dic, CLB23_dic$Core != 13.1), aes(x = DIC_original_mM, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "Depth (cm)", x = "DIC (mM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_porosity_plot <- ggplot(data = CLB23_porosity, aes(x = Porosity, y = Depth)) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(y = "Depth (cm)", x = "Porosity", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_hydrogen.summary <- CLB23_hydrogen %>% group_by(Depth, Core) %>%
  summarize(xmin = (H2_aq_nM - stdev),
            xmax = (H2_aq_nM + stdev),
            xmean = H2_aq_nM) %>%
  mutate(xmin = case_when(xmin < 0 ~ 0,
                          xmin > 0 ~ xmin))

CLB23_hydrogen_plot <- ggplot(data = CLB23_hydrogen.summary, aes(x = xmean, y = Depth)) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), color = "#77797a") +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(y = "Depth (cm)", x = "Hydrogen (nM)", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8)) +
  xlim(0,3.6)

CLB23_cellcounts.summary <- CLB23_cellcounts %>% group_by(Depth, Core) %>%
  summarize(xmin = (avg_cells_ml - std_dev_cells_ml),
            xmax = (avg_cells_ml + std_dev_cells_ml),
            xmean = avg_cells_ml) %>%
  mutate(xmin = case_when(xmin < 0 ~ 0,
                          xmin > 0 ~ xmin))

CLB23_cellcounts_plot <- ggplot(data = CLB23_cellcounts.summary, aes(x = xmean, y = Depth)) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), color = "#77797a") +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_x_log10(labels=trans_format("log10", math_format(10^.x))) +
  labs(y = "Depth (cm)", x = "Cells/g", shape = "Core", color = "Core") +
  theme(axis.title=element_text(size=8))

CLB23_geochemistry_main <- ggarrange(CLB23_methane_plot, CLB23_hydrogen_plot + labs(y = NULL), CLB23_d13CH4_plot + labs(y = NULL), CLB23_d13CO2_plot + labs(y = NULL), CLB23_dic_plot, CLB23_sulfate_plot + labs(y = NULL), CLB23_sulfide_plot + labs(y = NULL), CLB23_cellcounts_plot + labs(y = NULL), nrow = 2, ncol = 4,
                                labels = c("A", "B", "C", "D", "E", "F", "G", "H"), hjust = c(-.8, .3, .3, .3, -.8, .3, .3, .3), common.legend = TRUE)
CLB23_geochemistry_supplemental <- ggarrange(CLB23_porosity_plot, nrow = 1, ncol = 1,
                                     labels = c("A"), hjust = -.8, common.legend = TRUE)
setwd("~/Desktop/CLB23/Figures")

tiff("CLB23_geochemistry_main.tiff", units="in", width=8, height=6.4, res=300)
CLB23_geochemistry_main
dev.off()
tiff("CLB23_geochemistry_supplemental.tiff", units="in", width=4, height=4, res=300)
CLB23_geochemistry_supplemental
dev.off()

setwd("~/Desktop/CLB23/Data")

#variance calculations

WOR19_to_15_cmbsf <- subset.data.frame(WOR19, Depth < 15)
var(WOR19_to_15_cmbsf$H2_nM)

WOR19_below_20_cmbsf <- subset.data.frame(WOR19, Depth > 20)
var(WOR19_below_20_cmbsf$H2_nM)

CLB23_to_30_cmbsf <- subset.data.frame(CLB23_hydrogen, Depth < 30)
var(CLB23_to_30_cmbsf$H2_aq_nM)

CLB23_to_15_cmbsf <- subset.data.frame(CLB23_hydrogen, Depth < 15)
var(CLB23_to_15_cmbsf$H2_aq_nM)

##############################
#       thermodynamics       #
##############################

#plotting thermodynamic plots for methane production and sulfate reduction

#calculating ΔG° value for hydrogenotrophic methanogenesis
basis(c("H2O", "HCO3-", "H2", "H+"), c(1, -1, -4, -1), c("liq", "aq", "aq", "aq"))
subcrt("CH4", 1, "aq", T = 24, P = 1) #sample conditions T = 24 °C and P  = 1 atm
#activity coefficents
subcrt("HCO3-", T = 24, IS = 0.7, P = 1)

#calculating ΔG° value for H2-based sulfate reduction
subcrt(c("HS-", "H2O", "H2", "SO4-2"), c(1, 4, -4, -1), c("aq", "liq", "aq", "aq"), T = 24, P = 1)
#activity coefficents
subcrt("HS-", T = 24, IS = 0.7, P = 1)
subcrt("HSO4-", T = 24, IS = 0.7, P = 1)

#calculating ΔG° value for acetate-based sulfate reduction
subcrt(c("acetate", "SO4-2", "HS-", "HCO3-"), c(-1, -1, 1, 2), c("aq", "aq", "aq", "aq"), T = 24, P = 1)
#activity coefficents
subcrt("CH3COO-",T = 24, IS = 0.7, P = 1)

#calculating ΔG° value for propionate-based sulfate reduction
subcrt(c("propanoate", "SO4-2", "H+", "HS-", "HCO3-"), c(-4, -7, 1, 7, 12), c("aq", "aq", "aq", "aq", "aq"), T = 24, P = 1)
#activity coefficents
subcrt("propanoate",T = 24, IS = 0.7, P = 1)

#calculating ΔG° value for butyrate-based sulfate reduction
subcrt(c("n-butanoate", "SO4-2", "H+", "HS-", "HCO3-"), c(-2, -5, 1, 5, 8), c("aq", "aq", "aq", "aq", "aq"), T = 24, P = 1)
#activity coefficents
subcrt("n-butanoate",T = 24, IS = 0.7, P = 1)

#combined deltaG plots

CLB23_deltaG_methanogenesis_plot <- ggplot(data = CLB23_deltaG_methanogenesis, aes(x = deltaG_kj_mol, y = Depth)) +
  geom_hline(yintercept = 30, linetype = "dotdash", color = "#7CAE00", alpha = .5) +
  geom_hline(yintercept = 40, linetype = "dotdash", color = "#C77CFF", alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(y = "Depth (cm)", x = bquote("Δ"~italic(G[r])* ~ "kJ" ~ "mol"^-1 ~ CH[4]), shape = "Core", color = "Core", title = "Methanogenesis") +
  theme(axis.title=element_text(size=8))

CLB23_deltaG_h2_sulfatereduction_plot <- ggplot(data = CLB23_deltaG_sulfatereduction, aes(x = ΔG_hydrogen, y = Depth)) +
  geom_hline(yintercept = 30, linetype = "dotdash", color = "#7CAE00", alpha = .5) +
  geom_hline(yintercept = 40, linetype = "dotdash", color = "#C77CFF", alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
  geom_point(aes(shape = as.character(Core), color = as.character(Core)), size = 2.5) +
  scale_y_continuous(trans = "reverse") +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(y = "Depth (cm)", x = bquote("Δ"~italic(G[r])* ~ "kJ" ~ "mol"^-1 ~ SO[4]^-2), shape = "Core", color = "Core", title = "Sulfate reduction") +
  theme(axis.title=element_text(size=8))


CLB23_deltaG_main <- ggarrange(CLB23_deltaG_methanogenesis_plot, 
                               CLB23_deltaG_h2_sulfatereduction_plot + labs(y = NULL), nrow = 1, ncol = 2,
                                         labels = c("A", "B"), hjust = c(-.8, .3), common.legend = TRUE)

setwd("~/Desktop/CLB23/Figures")
tiff("CLB23_deltaG_main_onlyH2.tiff", units="in", width=4, height=4, res=300)
CLB23_deltaG_main
dev.off()
setwd("~/Desktop/CLB23/Data")

##############################
#            WOR             #
##############################

#creating WOR summary - methane, hydrogen, thermodynamics

WOR19_methane_plot <- ggplot(data = WOR19, aes(x = ch4_real, y = Depth)) +
  geom_point(aes(), size = 1) +
  scale_y_continuous(trans = "reverse") +
  labs(y = "Depth (cm)", x = "Methane (mM)", title = "") +
  theme(axis.title=element_text(size=8))

WOR19_hydrogen_plot <- ggplot(data = WOR19, aes(x = H2_nM, y = Depth)) +
  geom_point(aes(), size = 1) +
  scale_y_continuous(trans = "reverse") +
  labs(y = "Depth (cm)", x = "Hydrogen (nM)", title = "") +
  theme(axis.title=element_text(size=8)) +
  xlim(0,3.6)

WOR19_deltaG_methanogenesis_plot <- ggplot(data = WOR19, aes(x = deltaG_kj_mol, y = Depth)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
  geom_point(aes(), size = 1) +
  scale_y_continuous(trans = "reverse") +
  labs(y = "Depth (cm)", x = bquote("Δ"~italic(G[r])* ~ "kJ" ~ "mol"^-1 ~ CH[4]), title = "Methanogenesis") +
  theme(axis.title=element_text(size=8))

WOR19_summary_main <- ggarrange(WOR19_methane_plot, 
                               WOR19_hydrogen_plot + labs(y = NULL),
                               WOR19_deltaG_methanogenesis_plot + labs(y = NULL), nrow = 1, ncol = 3,
                               labels = c("A", "B", "C"), hjust = c(-.8, .3, .3), common.legend = TRUE)

setwd("~/Desktop/CLB23/Figures")
tiff("WOR19_summary_main.tiff", units="in", width=8, height=3.2, res=300)
WOR19_summary_main
dev.off()
setwd("~/Desktop/CLB23/Data")

#H2 leaktest for methods - for supplementary materials

H2_leaktest_plot <- ggplot(data = H2_leaktest, aes(x = time_elapsed_hr, y = avg)) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#77797a") +
  geom_point(aes(), size = 2.5) +
  labs(y = "Average peak area (Hydrogen)", x = "Hours elapsed") +
  theme(axis.title=element_text(size=8))

setwd("~/Desktop/CLB23/Figures")
tiff("H2_leaktest.tiff", units="in", width=4, height=4, res=300)
H2_leaktest_plot
dev.off()
setwd("~/Desktop/CLB23/Data")
