library(fdapace)

# Set the number of subjects (N) and the
# number of measurements per subjects (M) 
N <- 400;
M <- 100;

# Define the time interval
s <- seq(0,1,length.out = M)

# Define the mean curve and 2 eigencomponents
meanFunct <- function(s) s + 10*exp(-(s-5)^2)

eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)

# Create FPC scores
Ksi <- matrix(rnorm(N*2), ncol=2);
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(c(5,2))
dim(Ksi)

# Create Y_true
yTrue <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

# Sparsify the data
ySparse <- Sparsify(yTrue, s, sparsity = c(1:10))
names(ySparse)
ySparse$Lt # time points
ySparse$Ly # observations

# Add noise
ySparse$yNoisy <- lapply(ySparse$Ly, function(x) x + 0.5*rnorm(length(x)))

# Fit the model
Sparse <- FPCA(ySparse$yNoisy, ySparse$Lt, optns=list(plot = TRUE, FVEthreshold=0.9))

Sparse$sigma2
Sparse$lambda
Sparse$phi
Sparse$xiEst
Sparse$bwMu
Sparse$bwCov

Sparse <- FPCA(ySparse$yNoisy, ySparse$Lt, optns=list(plot = TRUE, FVEthreshold=0.95))

Sparse <- FPCA(ySparse$yNoisy, ySparse$Lt, optns=list(plot = TRUE, FVEthreshold=0.99,methodMuCovEst="smooth"))

