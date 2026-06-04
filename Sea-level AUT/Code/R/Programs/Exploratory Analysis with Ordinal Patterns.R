library(StatOrdPattHxC)
library(ggplot2)
library(ggthemes)

## For emb=3,4,5,6

attach(sea_level_multistation)

### Embedding =3

OrdinalPatterns <- OPseq(PUYT[!is.na(PUYT)], emb = 3)
OrdinalPatternsProbabilities <- OPprob(PUYT[!is.na(PUYT)], emb = 3)
HC_PUYT_3 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(AUCT[!is.na(AUCT)], emb = 3)
OrdinalPatternsProbabilities <- OPprob(AUCT[!is.na(AUCT)], emb = 3)
HC_AUCT_3 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(RFRT[!is.na(RFRT)], emb = 3)
OrdinalPatternsProbabilities <- OPprob(RFRT[!is.na(RFRT)], emb = 3)
HC_RFRT_3 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(NAPT[!is.na(NAPT)], emb = 3)
OrdinalPatternsProbabilities <- OPprob(NAPT[!is.na(NAPT)], emb = 3)
HC_NAPT_3 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

dim3 <- rbind(HC_AUCT_3, HC_NAPT_3, HC_RFRT_3, HC_PUYT_3)
dim3 <- data.frame(dim3)
names(dim3) <- c("H", "C")

ggplot(dim3, aes(x=H, y=C)) +
  geom_point() +
  geom_line(data=subset(LinfLsup, Dimension==3 & Side=="Lower"), aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==3 & Side=="Upper"), aes(x=H, y=C)) +
  theme_tufte()

### Embedding = 6

OrdinalPatterns <- OPseq(AUCT[!is.na(AUCT)], emb = 6)
OrdinalPatternsProbabilities <- OPprob(AUCT[!is.na(AUCT)], emb = 6)
HC_AUCT_6 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(PUYT[!is.na(PUYT)], emb = 6)
OrdinalPatternsProbabilities <- OPprob(PUYT[!is.na(PUYT)], emb = 6)
HC_PUYT_6 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(RFRT[!is.na(RFRT)], emb = 6)
OrdinalPatternsProbabilities <- OPprob(RFRT[!is.na(RFRT)], emb = 6)
HC_RFRT_6 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

OrdinalPatterns <- OPseq(NAPT[!is.na(NAPT)], emb = 6)
OrdinalPatternsProbabilities <- OPprob(NAPT[!is.na(NAPT)], emb = 6)
HC_NAPT_6 <- c(HShannon(OrdinalPatternsProbabilities), StatComplexity(OrdinalPatternsProbabilities))

dim6 <- rbind(HC_AUCT_6, HC_NAPT_6, HC_RFRT_6, HC_PUYT_6)
dim6 <- data.frame(dim6)
names(dim6) <- c("H", "C")

ggplot(dim6, aes(x=H, y=C)) +
  geom_point() +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Lower"), aes(x=H, y=C)) +
  geom_line(data=subset(LinfLsup, Dimension==6 & Side=="Upper"), aes(x=H, y=C)) +
  theme_tufte()

