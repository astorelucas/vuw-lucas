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

