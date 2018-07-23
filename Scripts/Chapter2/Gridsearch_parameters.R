library(BASiCS)
library(cowplot)

# Read in MCMCs run with different parameter combinations
MCMC.all <- readRDS("Google Drive File Stream/My Drive/BASiCS_add-on/Analysis/Results/Testing/Tdist/CD4_sample_MCMC_gridsearch.rds")

# Plot the fits
p.4_1.2_5 <- BASiCS_showFit(MCMC.all$CD4_4_1.2_5)
p.12_1.2_5 <- BASiCS_showFit(MCMC.all$CD4_12_1.2_5)
p.20_1.2_5 <- BASiCS_showFit(MCMC.all$CD4_20_1.2_5)
p.12_0.4_5 <- BASiCS_showFit(MCMC.all$CD4_12_0.4_5, variance = 0.4)
p.12_2_5 <- BASiCS_showFit(MCMC.all$CD4_12_2_5, variance = 2)
p.12_1.2_1 <- BASiCS_showFit(MCMC.all$CD4_12_1.2_1)
p.12_1.2_9 <- BASiCS_showFit(MCMC.all$CD4_12_1.2_9)

final <- plot_grid(p.4_1.2_5, p.12_1.2_5, p.20_1.2_5,
                   p.12_0.4_5, p.12_1.2_5, p.12_2_5,
                   p.12_1.2_1, p.12_1.2_5, p.12_1.2_9, ncol = 3)

ggsave("GitHub/Thesis/Chapter2/Figures/Fig_3.pdf", final, width = 12, height = 12)
