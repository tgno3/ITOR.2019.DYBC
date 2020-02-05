#########################################################################################################################
### Project  : A revised dynamic DEA model with budget constraints
### Script   : Dynamic DEA with budget constraints on GIT.R
### Contents : Operational assessment on fire departments
#########################################################################################################################

#########################################################################################################################
### Setting up Environment
#########################################################################################################################

# Load library
library(DJL)

# Load data and parameters
load("DYBC.RData")


#########################################################################################################################
### Numerical Example
#########################################################################################################################

# Table 1. Experimental data
table.1 <- cbind(df.e.io[, 3, 1] + df.e.io[, 3, 2] + df.e.Z.T, 
                 df.e.io[,  , 1],  df.e.io[, 3, 2] + df.e.Z.T, 
                 df.e.io[,  , 2],  df.e.Z.T)
print(table.1[, c(1:2, 4, 3, 5:6, 8, 7, 9)])


# Run
res.e.et <- dm.dynamic.et(df.e.io[,1,], df.e.io[,2,], df.e.io[,3,], df.e.Z.T)
data.frame(Sys.Eff = res.e.et$eff.s, T.Eff = res.e.et$eff.t, Lambda = res.e.et$lambda)


#########################################################################################################################
### Fire Department
#########################################################################################################################

# Table 2. Descriptive statistics
df.f.agg <- cbind(Total.Budget = rep(df.f.bg$Z.0 * 10^-6, 5),
                  df.f.2d[, c(id.x) + 1],
                  df.f.2d[, c(id.z, id.y[1]) + 1] * 10^-6,
                  df.f.2d[, c(id.y[2]) + 1, drop = F])
table.2 <- sapply(df.f.agg, function(x) c(Min  = min(x), 
                                          Med  = median(x), 
                                          Mean = mean(x), 
                                          Max  = max(x), 
                                          Std  = sd(x)))
noquote(format(round(t(table.2), 2), big.mark = ","))


# Table 3. Comparative results of efficiency
res.f.et <- dm.dynamic.et(df.f.3d[, id.x, ], df.f.3d[, id.y, ], df.f.3d[, id.z, ], df.f.bg$Z.T, rts)
res.f.bc <- dm.dynamic.bc(df.f.3d[, id.x, ], df.f.3d[, id.y, ], df.f.3d[, id.z, ], df.f.bg$Z.0, rts)
table.3  <- matrix(c(res.f.et$eff.s, res.f.bc$eff.s, res.f.et$eff.t, res.f.bc$eff.t), nrow(df.f.bg), 
                  dimnames = list(df.f.bg$DMU, 
                                  c("ET.Eff.sys", "LL.Eff.sys", paste0("ET.", 2012:2016), paste0("LL.", 2012:2016))))
#id.show  <- c(which(round(table.3[, 1], 6) == 1), which(rownames(table.3) == "Yeoju"))
print(round(table.3[, c(1, 2, 3, 8, 4, 9, 5, 10, 6, 11, 7, 12)], 4))


# How many system efficient DMUs?
apply(table.3[, 1:2], 2, function(x) sum(round(x, 8) == 1))


# Footnote 5
summary(lm((table.3[,1] - table.3[,2]) ~ df.f.bg$Z.T))


# Table 4. Namyangju
id.nyj  <- which(df.f.2d$DMU == "Namyangju")
table.4 <- rbind(df.f.agg[id.nyj, -1],
                 aggregate(df.f.agg[-id.nyj, -1], list(df.f.2d$Year[-id.nyj]), "mean")[, -1])
rownames(table.4) <- c(2012:2016, paste0("Avg.", 2012:2016))
round(t(table.4[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10),]), 2)

