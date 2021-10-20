function (results.table, plot.title = "Global shift in 3'UTR length", 
          do.ranksum.test = TRUE, return.plot = TRUE, do.plot = FALSE) 
{
  locations.res.table.up <- subset(results.table, FC_direction == 
                                     "Up")
  pos.upreg <- apply(as.matrix(locations.res.table.up[, c("SiteLocation", 
                                                          "NumSites")]), 1, function(x) {
                                                            relative_location(x[1], x[2])
                                                          })
  locations.res.table.down <- subset(results.table, FC_direction == 
                                       "Down")
  pos.downreg <- apply(as.matrix(locations.res.table.down[, 
                                                          c("SiteLocation", "NumSites")]), 1, function(x) {
                                                            relative_location(x[1], x[2])
                                                          })
  if (do.ranksum.test) {
    this.test <- wilcox.test(pos.upreg, pos.downreg)
    print("Wilcoxon Rank-sum test comparing relative peak locations for up- vs down-regulated peaks:")
    print(paste0("P-value = ", this.test$p.value))
    mytitle <- paste0("Wilcoxon Rank-sum test comparing \n", 
                      "relative peak locations for up- vs down-regulated peaks: \n", 
                      "P-value = ", this.test$p.value)
  }
  ggData <- data.frame(Peak_location = c(pos.upreg, pos.downreg), 
                       FC_direction = c(rep("Up", length(pos.upreg)), rep("Down", 
                                                                          length(pos.downreg))))
  ggData$FC_direction <- factor(ggData$FC_direction, levels = c("Up", 
                                                                "Down"))
  pl.density <- ggplot(ggData, aes(Peak_location, stat(count), 
                                   fill = FC_direction)) + geom_density(alpha = 0.8) + ylab("") + 
    theme_void(base_size = 18) + theme(axis.text = element_blank(), 
                                       axis.title = element_blank(), axis.ticks = element_blank()) + 
    ggtitle(mytitle) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    scale_fill_brewer(palette = "Set1")
  pl.histogram <- ggplot(ggData, aes(Peak_location, fill = FC_direction)) + 
    geom_histogram(position = position_dodge(), colour = "black", 
                   binwidth = 0.1, alpha = 0.8) + theme_classic(base_size = 18) + 
    xlab("Relative peak location") + ylab("Peak count") + 
    scale_fill_brewer(palette = "Set1") + theme(legend.position = "right") + 
    guides(fill = guide_legend(title = "Fold-change\ndirection", 
                               title.position = "top"))
  pl.combined <- cowplot::plot_grid(pl.density, pl.histogram, 
                                    ncol = 1, rel_heights = c(0.3, 0.8), axis = "lr", align = "v")
  if (do.plot) {
    plot(pl.combined)
  }
  if (return.plot) {
    return(pl.combined)
  }
}
