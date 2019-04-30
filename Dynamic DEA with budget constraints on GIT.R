#########################################################################################################################
### Project  : A revised dynamic DEA model with budget constraints
### Script   : Dynamic DEA with budget constraints on GIT.R
### Contents : Operational assessment on fire departments
#########################################################################################################################

#########################################################################################################################
### Setting up Environment
#########################################################################################################################

# Load library and functions
library(DJL)
source("dm.dynamic.et.R")

# Load data
df.f.bg <- read.csv(url("http://bit.ly/Budget4Fire"), header = T)
df.f.2d <- read.csv(url("http://bit.ly/Fire4Data"),   header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1, 11)], df.f.2d$Year, as.matrix))

# Parameter
id.t <- c(1)
id.x <- c(2:4)
id.y <- c(5:6)
id.z <- c(7)
rts  <- "vrs"


#########################################################################################################################
### Numerical Example
#########################################################################################################################

# Experimental data
df.e.io  <- array(c(2, 4, 8, 4, 1, 2, 2, 2, 3, 6, 12, 6,
                  5, 4, 3, 8, 1, 1, 1, 1, 5, 4,  3, 8),
                c(4, 3, 2), 
                dimnames = list(LETTERS[1:4], c("X", "Y", "z"), c("t1", "t2")))
df.e.Z.T <- array(c(1, 2, 3, 10), c(4, 1), dimnames = list(LETTERS[1:4], c("Z^T")))

# Run
res.e.et <- dm.dynamic.et(df.e.io[,1,], df.e.io[,2,], df.e.io[,3,], df.e.Z.T)
data.frame(Sys.Eff = res.e.et$eff.s, T.Eff = res.e.et$eff.t, Lambda = res.e.et$lambda)


#########################################################################################################################
### Fire Department
#########################################################################################################################

# Table 1. Descriptive statistics
df.f.agg <- cbind(Total.Budget = rep(df.f.bg$Z.0 * 10^-6, 5),
                  df.f.2d[, c(id.x) + 1],
                  df.f.2d[, c(id.z, id.y[1]) + 1] * 10^-6,
                  df.f.2d[, c(id.y[2]) + 1, drop = F])
table.1 <- sapply(df.f.agg, function(x) c(Min  = min(x), 
                                          Med  = median(x), 
                                          Mean = mean(x), 
                                          Max  = max(x), 
                                          Std  = sd(x)))
noquote(format(round(t(table.1), 2), big.mark = ","))

# Table 2. Comparative results of efficiency
res.f.et <- dm.dynamic.et(df.f.3d[, id.x, ], df.f.3d[, id.y, ], df.f.3d[, id.z, ], df.f.bg$Z.T, rts)
res.f.bc <- dm.dynamic.bc(df.f.3d[, id.x, ], df.f.3d[, id.y, ], df.f.3d[, id.z, ], df.f.bg$Z.0, rts)
table.2  <- matrix(c(res.f.et$eff.s, res.f.bc$eff.s, res.f.et$eff.t, res.f.bc$eff.t), nrow(df.f.bg), 
                  dimnames = list(df.f.bg$DMU, 
                                  c("ET.Eff.sys", "LL.Eff.sys", paste0("ET.", 2012:2016), paste0("LL.", 2012:2016))))
print(round(table.2[, c(1, 2, 3, 8, 4, 9, 5, 10, 6, 11, 7, 12)], 4))

# How many system efficient DMUs?
apply(table.2[, 1:2], 2, function(x) sum(round(x, 8) == 1))

# Footnote 5
summary(lm((table.2[,1] - table.2[,2]) ~ df.f.bg$Z.T))

# Table 3. Namyangju
id.nyj  <- which(rownames(table.2) == "Namyangju")
table.3 <- rbind(df.f.agg[id.nyj + 33 * 0:4, -1],
                 aggregate(df.f.agg[, -1], list(df.f.2d$Year), "mean")[, -1])
rownames(table.3) <- c(2012:2016, paste0("Avg.", 2012:2016))
round(t(table.3[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10),]), 2)

