#Figure 2F
rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ggrepel)

data_dir <- "/data/zdh/ZDH_Spatial/Kidney_Spatial_Project/ceshi/Figure/20260922"
csv_path <- file.path(data_dir, "PF4_Mono_vs_rest_DEGs_clean.csv")
deg_raw <- read.csv(csv_path)

deg <- deg_raw %>%
  mutate(
    p_val_adj_clean = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_clean)
  ) %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(rank = row_number())

palette_colors <- rev(brewer.pal(n = 11, name = "RdYlBu"))
col_palette <- colorRampPalette(palette_colors)(100)

top_pos_genes <- c("GP9", "CMTM5", "ITGB3", "TUBB1", "PPBP", "ITGA2B", "PF4", "CXCL8")
top_pos <- deg %>% filter(Gene %in% top_pos_genes)

top_neg_genes <- c("IL18R1", "STARD13", "KCNMA1", "IL1R2", "VSIG4")
top_neg <- deg %>% filter(Gene %in% top_neg_genes)

max_rank <- max(deg$rank)

p <- ggplot(deg, aes(x = rank, y = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.5) +
  geom_point(aes(color = neg_log10_padj, size = neg_log10_padj), alpha = 0.85) +
  scale_size_continuous(range = c(1.8, 3.8), guide = "none") +
  scale_color_gradientn(
    colors = col_palette,
    guide = guide_colorbar(
      title = expression(-log[10](italic(P)[adj])),
      title.position = "top",
      barwidth = unit(0.7, "lines"),
      barheight = unit(5.5, "lines"),
      ticks.colour = "black"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::comma) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Ranked DEGs",
    y = expression(bold(avg_log[2] ~ "Fold Change")),
    title = expression(bold(PF4^"+" * PPBP^"+" ~ "Mono vs. Rest Monocytes"))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", linewidth = 1.1, fill = NA),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = "right"
  )

p_final <- p +
  geom_text_repel(
    data = top_pos,
    aes(label = Gene),
    color = "#4A6B82",
    fontface = "bold",
    size = 4.0,
    direction = "y",
    hjust = 0,
    nudge_x = max_rank * 0.12,
    segment.size = 0.35,
    segment.color = "grey60",
    min.segment.length = 0,
    box.padding = 0.35,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = top_neg,
    aes(label = Gene),
    color = "#4A6B82",
    fontface = "bold",
    size = 4.0,
    direction = "y",
    hjust = 1,
    nudge_x = -max_rank * 0.10,
    segment.size = 0.35,
    segment.color = "grey60",
    min.segment.length = 0,
    box.padding = 0.4,
    force = 2.5,
    max.overlaps = Inf
  )

print(p_final)

pdf_out <- file.path(data_dir, "Figure_PF4_Mono_Raw_Truth.pdf")
ggsave(pdf_out, plot = p_final, width = 5.2, height = 7.5, dpi = 300)