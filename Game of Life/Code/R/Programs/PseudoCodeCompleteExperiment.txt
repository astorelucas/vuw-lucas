# Structure

Input data:
  Seed
  k: the image dimension is 2^k (assume 256, then generalise)
  T: number of epochs
  D: embedding dimension(s)
  
  
Intialize map0
For t in 1:T
  map <- map from map t-1
  // point in HxC from map
  time.series <- NULL
  For i in 1:253:4
    For j in 1:253:4
      window.binary <- as.vector(map[i:(i+3), j:(j+3)])
      window.real <- binary.to.real(window.binary)
      time.series <- c(time.series, window.real)
    Endfor
  Endfor
  hc <- HxC(time.series)
Endfor