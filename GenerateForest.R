# 1. Install & load
library(dplyr)
library(forestplot)
library(grid)    # for unit() and gpar()

# 2. Read only the 3rd,5th,6th,7th cols from TSV
dat <- read.delim("PanCancer_hallmarkgenes_OS_survival_results.tsv", sep="\t", 
                  stringsAsFactors=FALSE, check.names=FALSE) %>%
  select(3,6,7,8) %>%
  setNames(c("Hallmark","HR","CI_low","CI_high"))

# 3. Ensure numeric
dat <- dat %>%
  mutate(
    HR       = as.numeric(HR),
    CI_low   = as.numeric(CI_low),
    CI_high  = as.numeric(CI_high)
  )

# 4. Aggregate by Hallmark
dat_sum <- dat %>%
  group_by(Hallmark) %>%
  summarize(
    mean_HR   = mean(HR,      na.rm=TRUE),
    mean_low  = mean(CI_low,  na.rm=TRUE),
    mean_high = mean(CI_high, na.rm=TRUE),
    .groups   = "drop"
  )

# 5. Build tabletext
tabletext <- rbind(
  c("Hallmark","Mean HR","95% CI"),
  cbind(
    dat_sum$Hallmark,
    sprintf("%.2f", dat_sum$mean_HR),
    paste0("[", sprintf("%.2f", dat_sum$mean_low), ", ",
           sprintf("%.2f", dat_sum$mean_high), "]")
  )
)

# 6. Prepare numeric vectors
mean_vec  <- c(NA, dat_sum$mean_HR)
lower_vec <- c(NA, dat_sum$mean_low)
upper_vec <- c(NA, dat_sum$mean_high)

# 7. Plot
forestplot(
  labeltext  = tabletext,
  mean       = mean_vec,
  lower      = lower_vec,
  upper      = upper_vec,
  zero       = 1,
  boxsize    = 0.4,
  lineheight = unit(8, "mm"),
  colgap     = unit(4, "mm"),
  hrzl_lines = list("1" = gpar(lwd=1)),
  txt_gp     = fpTxtGp(
    label = gpar(cex=0.8),
    ticks = gpar(cex=0.7),
    xlab  = gpar(cex=0.9)
  ),
  xlab       = "Hazard Ratio (HR)",
  title      = "Average HR by Hallmark Type",
  lwd.ci     = 4   # ← increase from default (usually 1) to e.g. 2
)


