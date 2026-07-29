install.packages(c(
  "EBImage",
  "igraph",
  "glcm",
  "moments",
  "pracma",
  "imager",
  "e1071"
))

library(EBImage)
library(igraph)
library(glcm)
library(pheatmap)

install.packages("igraph", type = "source")
install.packages("BiocManager")
BiocManager::install("EBImage")

density_feature <- function(img){
  
  mean(img)
  
}

neighbor_histogram <- function(img){
  
  nr <- nrow(img)
  nc <- ncol(img)
  
  counts <- numeric(9)
  
  for(i in 2:(nr-1)){
    for(j in 2:(nc-1)){
      
      n <- sum(img[(i-1):(i+1),
                   (j-1):(j+1)]) -
        img[i,j]
      
      counts[n+1] <- counts[n+1] + 1
    }
  }
  
  counts/sum(counts)
  
}

connected_features <- function(img){
  
  cc <- bwlabel(Image(img))
  
  areas <- table(cc)
  
  areas <- areas[names(areas)!="0"]
  
  if(length(areas)==0){
    
    return(c(
      components=0,
      largest=0,
      mean_area=0
    ))
    
  }
  
  c(
    components=length(areas),
    largest=max(areas),
    mean_area=mean(areas)
  )
  
}

euler_feature <- function(img){
  
  cc <- bwlabel(Image(img))
  
  components <- length(unique(cc))-1
  
  components
  
}

glcm_features <- function(img){
  
  contrast <- glcm(
    img,
    window=c(3,3),
    statistics="contrast"
  )
  
  homogeneity <- glcm(
    img,
    window=c(3,3),
    statistics="homogeneity"
  )
  
  
  correlation <- glcm(
    img,
    window=c(3,3),
    statistics="correlation"
  )
  
  
  c(
    contrast = mean(contrast, na.rm=TRUE),
    homogeneity = mean(homogeneity, na.rm=TRUE),
    correlation = mean(correlation, na.rm=TRUE)
  )
  
}

local_density <- function(img){
  
  kernel <- matrix(1,3,3)
  kernel[2,2] <- 0
  
  convolve2d <- function(x,k){
    as.matrix(
      convolve(
        x,
        rev(k),
        type="filter"
      )
    )
  }
  
  convolve2d(img,kernel)
}

fft_features <- function(img){
  
  F <- fft(img)
  
  power <- Mod(F)^2
  
  c(
    
    fft_mean=mean(power),
    
    fft_sd=sd(as.vector(power)),
    
    fft_entropy=-sum(
      (power/sum(power))*
        log(power/sum(power)+1e-12)
    )
    
  )
  
}

boxcount <- function(img){
  
  sizes <- c(2,4,8,16,32)
  
  N <- numeric(length(sizes))
  
  for(k in seq_along(sizes)){
    
    s <- sizes[k]
    
    count <- 0
    
    for(i in seq(1,nrow(img),by=s))
      for(j in seq(1,ncol(img),by=s)){
        
        block <- img[
          i:min(i+s-1,nrow(img)),
          j:min(j+s-1,ncol(img))
        ]
        
        if(any(block==1))
          count <- count+1
      }
    
    N[k] <- count
  }
  
  fit <- lm(log(N)~log(1/sizes))
  
  coef(fit)[2]
  
}

spatial_entropy <- function(img){
  
  p <- table(img)/length(img)
  
  -sum(p*log2(p))
  
}

extract_features <- function(img){
  
  neigh <- neighbor_histogram(img)
  
  names(neigh) <- paste0("n",0:8)
  
  c(
    
    density=density_feature(img),
    
    connected_features(img),
    
    euler=euler_feature(img),
    
    spatial_entropy=spatial_entropy(img),
    
    fft_features(img),
    
    fractal_dimension=boxcount(img),
    
    glcm_features(img),
    
    neigh
    
  )
  
}

## Running:

epochs <- c(10)
n_replicates <- 50
features_list <- list()
results_list <- list()
i<-2
for (e in epochs){
  print(e)
  # Peguei o estado inicial e o resultado daquela que quero analisar
  features_list[[i]] <- ca_history10k_0.8[[1]]
  
  results <- paste0("Data/results/", e ,"K/results_", e, "k_0.8_emb5_16bits.csv")
  r <- read.csv(results)
  results_list[[i]] <- r
  
  for (n in 2:n_replicates){
    ca_simulation2 <- paste0("Data/ca-simulations/", e ,"K/0.8_100execucoes/ca_history_", e, "k_0.8_", n, ".rds")
    results <- paste0("Data/results/", e ,"K/emb5_0.8_16bits/results_", e, "k_0.8_emb5_16bits_", n, ".csv")
    
    img <- readRDS(ca_simulation2)
    r2 <- read.csv(results)
    features_list[[i]] <- img[[1]]
    results_list[[i]] <- r2
    
    print(n)
    i<-i+1
  }
}

features_list[[1]] <- ca_history10k_0.8[[1]]

results <- paste0("Data/results/10K/results_10k_0.8_emb5_16bits.csv")
r <- read.csv(results)
results_list[[1]] <- r

feature_matrix <- do.call(
  rbind,
  lapply(features_list, extract_features)
)

## Calcular média final C, média final H

feature_matrix
features_scaled <- scale(feature_matrix)
dist_matrix <- dist(
  features_scaled,
  method="euclidean"
)

pheatmap(
  as.matrix(dist_matrix),
  main="Distance between initial configurations"
)

features_clean <- features_scaled[, colSums(is.na(features_scaled)) == 0]

features_filtered <- features_clean[
  ,
  apply(features_clean,2,sd) > 0
]

pca <- prcomp(
  features_filtered,
  center=TRUE,
  scale.=TRUE
)

rownames(features_filtered) <- paste0("", 1:nrow(features_filtered))

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  
  configuration = rownames(features_filtered)
)

pca_df

## Plotting the PCA + Final C + Convergence_time

ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    #color=final_C,
    #size=convergence_time
  )
)+
  geom_point(
    alpha=0.8
  )+
  scale_color_viridis_c(
    name="Final C"
  )+
  scale_size_continuous(
    name="Convergence time"
  )+
  theme_bw(base_size=15)+
  labs(
    title="Initial Configuration Space and Emergent Complexity",
    subtitle="PCA projection of spatial descriptors",
    x="PC1",
    y="PC2"
  )
