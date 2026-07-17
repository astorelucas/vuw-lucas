library(readxl)
library(StatOrdPattHxC)
library(ggplot2)
library(ggthemes)
library(ggrepel)

sea_level_multistation <- read_csv("Sea-level AUT/Dataset/sea_level_multistation.csv", show_col_types = FALSE)

## For emb = 3,4,6

attach(sea_level_multistation)
stations <- c("AUCT", "CHIT", "CHST", "CPIT", "GBIT", "GIST", "JBTT", "KAIT", 
              "LOTT", "MNKT", "NAPT", "NBRT", "NCPT", "OTAT", "PUYT", "RBCT", 
              "RFRT", "TAUT", "TPTT", "WLGT")

### Embedding = 3

results_list <- list()

for (stn in stations) {
  x <- sea_level_multistation[[stn]]
  x <- x[!is.na(x)]
  
  OrdinalPatterns <- OPseq(x, emb = 3)
  OrdinalPatternsProbabilities <- OPprob(x, emb = 3)
  
  results_list[[stn]] <- c(HShannon(OrdinalPatternsProbabilities), 
                           StatComplexity(OrdinalPatternsProbabilities))
}

dim3 <- do.call(rbind, results_list)
dim3 <- data.frame(dim3)
names(dim3) <- c("H", "C")
rownames(dim3) <- stations

ggplot(dim3, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==3 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==3 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3,
                  fontface="bold",
                  max.overlaps=Inf,
                  box.padding=0.4,
                  point.padding=0.3,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(0.0, 1)) +
  scale_y_continuous(limits=c(0.0, 0.35)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 3",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))


### Embedding = 4

results_list4 <- list()

for (stn in stations) {
  x <- sea_level_multistation[[stn]]
  x <- x[!is.na(x)]
  
  OrdinalPatterns <- OPseq(x, emb = 4)
  OrdinalPatternsProbabilities <- OPprob(x, emb = 4)
  
  results_list4[[stn]] <- c(HShannon(OrdinalPatternsProbabilities), 
                            StatComplexity(OrdinalPatternsProbabilities))
}

dim4 <- do.call(rbind, results_list4)
dim4 <- data.frame(dim4)
names(dim4) <- c("H", "C")
rownames(dim4) <- stations

ggplot(dim4, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==4 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==4 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3,
                  fontface="bold",
                  max.overlaps=Inf,
                  box.padding=0.4,
                  point.padding=0.3,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(0.0, 1)) +
  scale_y_continuous(limits=c(0.0, 0.35)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 4",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

### Embedding = 5

results_list5 <- list()

for (stn in stations) {
  print(stn)
  x <- sea_level_multistation[[stn]]
  x <- x[!is.na(x)]
  
  OrdinalPatternsProbabilities <- OPprob(x, emb = 5)
  OrdinalPatternsProbabilities <- OrdinalPatternsProbabilities / sum(OrdinalPatternsProbabilities)  
  h <- HShannon(OrdinalPatternsProbabilities)
  c <- StatComplexity(OrdinalPatternsProbabilities)
  results_list5[[stn]] <- c(h, c)
}


dim5 <- do.call(rbind, results_list5)
dim5 <- data.frame(dim5)
names(dim5) <- c("H", "C")
rownames(dim5) <- stations

ggplot(dim5, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==5 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==5 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3,
                  fontface="bold",
                  max.overlaps=Inf,
                  box.padding=0.4,
                  point.padding=0.3,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(min(dim5$H), max(dim5$H))) +
  scale_y_continuous(limits=c(0.3, 0.4)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 5",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))


### Embedding = 6

results_list6 <- list()

for (stn in stations) {
  print(stn)
  x <- sea_level_multistation[[stn]]
  x <- x[!is.na(x)]
  
  OrdinalPatternsProbabilities <- OPprob(x, emb = 6)
  OrdinalPatternsProbabilities <- OrdinalPatternsProbabilities / sum(OrdinalPatternsProbabilities)  
  h <- HShannon(OrdinalPatternsProbabilities)
  c <- StatComplexity(OrdinalPatternsProbabilities)
  results_list6[[stn]] <- c(h, c)
}


dim6 <- do.call(rbind, results_list6)
dim6 <- data.frame(dim6)
names(dim6) <- c("H", "C")
rownames(dim6) <- stations

ggplot(dim6, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(0., 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 6",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

dim6$label <- ifelse(dim6$stations %in% c("PUYT", "OTAT"), dim6$stations, NA)

dim6$label <- ifelse(dim6$stations %in% c("PUYT", "OTAT"), dim6$stations, NA)

ggplot(dim6, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=label, color=stations),
                  size=3,
                  fontface="bold",
                  na.rm=TRUE,
                  max.overlaps=Inf,
                  box.padding=2,
                  point.padding=0.9,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 6",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

# INSERT PLOT

library(patchwork)

# Main plot (your existing code, saved to a variable)
main_plot <- ggplot(dim6, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=label, color=stations),
                  size=3,
                  fontface="bold",
                  na.rm=TRUE,
                  max.overlaps=Inf,
                  box.padding=2,
                  point.padding=0.9,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 6",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

station_levels <- sort(unique(dim6$stations))
station_colors <- setNames(
  colorRampPalette(c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(length(station_levels)),
  station_levels
)

dim6_zoom <- dim6 %>% 
  filter(stations %in% c("PUYT", "OTAT")) %>% 
  droplevels()
# Subset just the two stations of interest
zoom_stations <- c("PUYT", "OTAT")
dim6_zoom <- subset(dim6, stations %in% zoom_stations)

# Padding around the zoomed points so they're not flush with the edges
pad_x <- diff(range(dim6_zoom$H)) * 0.3
pad_y <- diff(range(dim6_zoom$C)) * 0.3
xlim_zoom <- range(dim6_zoom$H) + c(-pad_x, pad_x)
ylim_zoom <- range(dim6_zoom$C) + c(-pad_y, pad_y)

# Inset plot: same style, cropped to the zoom region
inset_plot <- ggplot(dim6_zoom, aes(x=H, y=C)) +
  geom_point(aes(color=stations), size=3, alpha=0.9) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3,
                  fontface="bold",
                  box.padding=0.5,
                  point.padding=0.3,
                  segment.color="grey60",
                  segment.size=0.4,
                  seed=42) +
  scale_color_manual(values=colorRampPalette(
    c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(20)) +
  coord_cartesian(xlim=xlim_zoom, ylim=ylim_zoom) +
  theme_tufte(base_size = 8) +
  theme(legend.position="none",
        axis.title=element_blank(),
        plot.background=element_rect(color="grey40", fill="white", linewidth=0.5))

# Combine: place inset in the top-right corner (adjust left/right/top/bottom as needed)
main_plot + inset_element(inset_plot, 
                          left = 0.62, right = 0.98, 
                          bottom = 0.62, top = 0.98)

dim6$stations <- stations   # attach it once, assuming same order/length



library(ggplot2)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(dplyr)

# --- Consistent color mapping across main + inset ---
station_levels <- sort(unique(dim6$stations))
station_colors <- setNames(
  colorRampPalette(c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(length(station_levels)),
  station_levels
)

# --- Label only PUYT and OTAT ---
dim6$label <- ifelse(dim6$stations %in% c("PUYT", "OTAT"), dim6$stations, NA)

# --- Zoom box coordinates (around PUYT/OTAT cluster, with small padding) ---
dim6_zoom <- dim6 %>% 
  filter(stations %in% c("PUYT", "OTAT")) %>% 
  droplevels()

pad_x <- diff(range(dim6_zoom$H)) * 0.3
pad_y <- diff(range(dim6_zoom$C)) * 0.3
# Guard against zero-width padding if H/C are near-identical
pad_x <- max(pad_x, 0.002)
pad_y <- max(pad_y, 0.0005)

xlim_zoom <- range(dim6_zoom$H) + c(-pad_x, pad_x)
ylim_zoom <- range(dim6_zoom$C) + c(-pad_y, pad_y)

# --- Main plot ---
main_plot <- ggplot(dim6, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  # Zoom box around PUYT/OTAT
  annotate("rect", xmin=xlim_zoom[1], xmax=xlim_zoom[2], 
           ymin=ylim_zoom[1], ymax=ylim_zoom[2],
           fill=NA, color="grey30", linewidth=0.5, linetype="dotted") +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=label, color=stations),
                  size=3,
                  fontface="bold",
                  na.rm=TRUE,
                  max.overlaps=Inf,
                  box.padding=1,
                  point.padding=0.5,
                  force=10,
                  force_pull=0.5,
                  min.segment.length=0,
                  segment.color="grey60",
                  segment.size=0.4,
                  nudge_x=0.05,
                  nudge_y=0.08,
                  seed=42) +
  scale_color_manual(values=station_colors) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(title="HxC Complexity-Entropy Plane",
       subtitle="Sea level stations — Embedding dimension 6",
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

# --- Inset plot: same color mapping as main plot ---
inset_plot <- ggplot(dim6_zoom, aes(x=H, y=C)) +
  geom_point(aes(color=stations), size=3, alpha=0.9) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3, fontface="bold",
                  box.padding=0.5, point.padding=0.3,
                  segment.color="grey60", segment.size=0.4, seed=42) +
  scale_color_manual(values=station_colors) +
  coord_cartesian(xlim=xlim_zoom, ylim=ylim_zoom) +
  theme_tufte(base_size = 8) +
  theme(legend.position="none",
        axis.title=element_blank(),
        plot.background=element_rect(color="grey40", fill="white", linewidth=0.5))

# --- Combine ---
main_plot + inset_element(inset_plot, 
                          left = 0.62, right = 0.98, 
                          bottom = 0.62, top = 0.98)








#---

library(ggplot2)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(dplyr)

# --- Consistent color mapping across main + inset ---
station_levels <- sort(unique(dim6$stations))
station_colors <- setNames(
  colorRampPalette(c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"))(length(station_levels)),
  station_levels
)

# --- Label only PUYT and OTAT ---
dim6$label <- ifelse(dim6$stations %in% c("PUYT", "OTAT"), dim6$stations, NA)

# --- Zoom region (actual data bounds around PUYT/OTAT) ---
dim6_zoom <- dim6 %>% 
  filter(stations %in% c("PUYT", "OTAT")) %>% 
  droplevels()

pad_x <- max(diff(range(dim6_zoom$H)) * 0.3, 0.002)
pad_y <- max(diff(range(dim6_zoom$C)) * 0.3, 0.0005)

xlim_zoom <- range(dim6_zoom$H) + c(-pad_x, pad_x)
ylim_zoom <- range(dim6_zoom$C) + c(-pad_y, pad_y)

# --- Box position on the main plot: shifted left of the actual cluster ---
# --- Box size: fixed, visible dimensions (independent of tiny zoom range) ---
box_size_x <- 0.10   # width of the visible box, in main-plot H units
box_size_y <- 0.10   # height of the visible box, in main-plot C units

# Center the box on the actual PUYT/OTAT cluster, then shift left
cluster_center_x <- mean(xlim_zoom)
cluster_center_y <- mean(ylim_zoom)

box_shift <- 0.8    # how far left of the cluster to draw the box

box_xmax <- cluster_center_x - box_shift + box_size_x/2
box_xmin <- box_xmax - box_size_x
box_ymin <- cluster_center_y - box_size_y/2
box_ymax <- cluster_center_y + box_size_y/2

# --- Main plot ---
main_plot <- ggplot(dim6, aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6, linetype="dashed",
            inherit.aes=FALSE) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), 
            aes(x=H, y=C), color="grey40", linewidth=0.6,
            inherit.aes=FALSE) +
  # Zoom box, shifted left of the actual points
  annotate("rect", xmin=box_xmin, xmax=box_xmax, ymin=box_ymin, ymax=box_ymax,
           fill=NA, color="grey30", linewidth=0.5, linetype="dotted") +
  # Connector segments: box corners -> inset location (top-right of plot area)
  annotate("segment", x=box_xmax, y=box_ymax, xend=0.62, yend=0.98,
           color="grey50", linewidth=0.4, linetype="dotted") +
  annotate("segment", x=box_xmax, y=box_ymin, xend=0.62, yend=0.62,
           color="grey50", linewidth=0.4, linetype="dotted") +
  geom_point(aes(color=stations), size=3, alpha=0.85) +
  geom_text_repel(aes(label=label, color=stations),
                  size=3,
                  fontface="bold",
                  na.rm=TRUE,
                  max.overlaps=Inf,
                  box.padding=0.8,
                  point.padding=0.5,
                  force=10,
                  force_pull=0.5,
                  min.segment.length=0,
                  segment.color="grey60",
                  segment.size=0.4,
                  nudge_x=0.05,
                  nudge_y=0.08,
                  seed=42) +
  scale_color_manual(values=station_colors) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(
       x="Permutation Entropy (H)",
       y="Statistical Complexity (C)") +
  theme_tufte() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40", size=10),
        axis.title=element_text(size=11))

# --- Inset plot: same color mapping as main plot ---
inset_plot <- ggplot(dim6_zoom, aes(x=H, y=C)) +
  geom_line(color="grey50", linewidth=0.5) +
  geom_point(aes(color=stations), size=3, alpha=0.9) +
  geom_text_repel(aes(label=stations, color=stations),
                  size=3, fontface="bold",
                  box.padding=0.5, point.padding=0.3,
                  segment.color="grey60", segment.size=0.4, seed=42) +
  scale_color_manual(values=station_colors) +
  coord_cartesian(xlim=xlim_zoom, ylim=ylim_zoom) +
  theme_tufte(base_size = 8) +
  theme(legend.position="none",
        axis.title=element_blank(),
        plot.background=element_rect(color="grey40", fill="white", linewidth=0.5))

# --- Combine ---
main_plot + inset_element(inset_plot, 
                          left = 0.05, right = 0.4, 
                          bottom = 0.52, top = 0.98)

