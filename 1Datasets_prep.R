packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", 
  "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
  "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
  "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
  "installr", "resemble", "prospectr", "magrittr", "doParallel", 
  "parallel", "foreach", "ggplot2", "tidyr", "pls")
# Install and load packages
# install.packages(packages_to_install, dependencies = TRUE)
# Load the installed packages
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}


# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1) {
  cal_statistics_result <- NULL
  for (i in seq(1, 40, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

#Stratified Function
stratified <- function(df, group, size, select = NULL, 
                       replace = FALSE, bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)
  
  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}

# library(readxl)
# library(writexl)
# #library(ggplot2)
# library(dplyr)
# library(reshape2)
# library(GGally)
# library(gridExtra)
# library(grid)
# library(readr)
# library(corrplot)
# library(plot.matrix)
# #install.packages("aricode")
# library(aricode)
# library(infotheo)
# library(heatmap3)
# library(pheatmap)
# library(lattice)
# #library(NPRED)- check appropriate version of R and install the same to continue
# library(NPRED)
# library(csvread)
# #install.packages("plotrix")
# library(soiltexture)
# library(stringr)
# library(installr)
# library(resemble)
# library(prospectr)
# library(magrittr)
# library(doParallel)
# library(parallel)
# library(foreach)
# library(ggplot2)
# library(pls)

##Loading the datasets##
data(NIRsoil)
dim(NIRsoil)
str(NIRsoil)
mode(NIRsoil)
attributes(NIRsoil)
mode(NIRsoil$spc)
attributes(NIRsoil$spc)

setwd("D:/Academics/PhD/SEM 11/Paper_3_Global2Local/Working_files/Refined/")
setwd("F:/CRT/Refined/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,9:length(foo)]))
attributes(foo1)
attributes(NIRsoil$spc)
dimnames(foo1)
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(400, 2500, by = 2))
attributes(foo1)
dimnames(foo1)
mode(foo1)
attributes(foo1$spc)
dimnames(foo1$spc)
mode(foo1$spc)
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Karnataka Data
foo = read_excel("Karnataka.xlsx", sheet = "Sheet1")
foo1 = as.data.frame(foo[,1:7])
foo2 = -log10(data.matrix(foo[,8:length(foo)]))
rownames(foo2) = as.character(seq(1, 497, by = 1))
colnames(foo2) = as.character(seq(400, 2500, by = 1))
#pre-process the spectra: resample it to a resolution of 2 nm 
old_wavs = foo2 %>% colnames() %>% as.numeric()
new_wavs = seq(400, 2500, by = 2)
foo1$spc = foo2 %>% 
  resample(wav = old_wavs, new.wav = new_wavs)
attributes(foo1)
dimnames(foo1)
mode(foo1)
attributes(foo1$spc)
dimnames(foo1$spc)
mode(foo1$spc)
kar = foo1
rm(foo,foo1,foo2)

#plot data
foo = OSSL
foo = kar
wavelengths = as.matrix(as.numeric(colnames(foo$spc)))
matplot(x = wavelengths, y = t(foo$spc), 
        xlab = "Wavelengths (nm)",
        ylab = "Absorbance",
        type = "l", lty = 1,lwd = 2, col = "#5177A133",
        ylim = c(0, 2),
        yaxt = "n",
        cex.main = 1.2,
        cex.lab = 1.2, 
        cex.axis = 1.2,
        font.axis = 2, 
        family = "sans",
        las = 1)
y_ticks <- seq(0, 2, by = 0.5)
axis(side = 2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = 1.2)

means = colMeans(foo$spc)
std_devs = apply(foo$spc, 2, sd)
percentiles = apply(foo$spc, 2, quantile, c(0, 0.25, 0.5, 0.75, 1))
spectra_data = data.frame(Wavelength = wavelengths, Mean = means, SD = std_devs)
ggplot(spectra_data, aes(x = Wavelength)) +
  geom_line(aes(y = Mean, color = "blue"), size = 1) +
  geom_line(aes(y = percentiles[5,], color = "red"), size = 1) +
  geom_line(aes(y = percentiles[1,], color = "green"), size = 1) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = "SD"), alpha = 0.5) +
  scale_color_manual(values = c("blue", "green", "red")) +
  scale_fill_manual(values = "lightblue", name = "1 SD") +
  labs(x = "Wavelength (nm)",y = "Absorbance") +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid = element_blank(),legend.position = "none",
        axis.ticks = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  coord_cartesian(ylim = c(0, 2))
rm(y_ticks, spectra_data, percentiles, means, std_devs)

#try later
foo_ = prcomp(foo$spec, scale. = TRUE)
foo_pca = as.data.frame(foo_$x[, 1:3])
# Scatter plot between PC1 and PC2
ggplot(spectral_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(x = "PC1", y = "PC2") +
  theme_minimal()

ggplot(spectral_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = stat(density)), alpha = 0.5) +  # Add scatter points with density-based color
  geom_density_2d(aes(fill = stat(density)), alpha = 0.5) +  # Add 2D density plot with density-based fill color
  scale_fill_viridis_c() +  # Use the Viridis color palette for fill color
  scale_color_viridis_c() +  # Use the Viridis color palette for scatter point color
  labs(x = "PC1", y = "PC2") +
  theme_minimal()

#try later ends

rm(foo)


#Data provision
IP = OSSL
IP1 = kar
#Component selection
property = "Clay"

X = IP$spc[!is.na(IP[,property]),]
Y = IP[,property][!is.na(IP[,property])]
IP1 = IP1[order(IP1[,property]), ]
IP1_training = IP1[c(TRUE, FALSE), ]
IP1_testing = IP1[c(FALSE, TRUE), ]
X1 = IP1_training$spc[!is.na(IP1_training[,property]),]
Y1 = IP1_training[,property][!is.na(IP1_training[,property])]
X2 = IP1_testing$spc[!is.na(IP1_testing[,property]),]
Y2 = IP1_testing[,property][!is.na(IP1_testing[,property])]

#First derivative calculation
X10 = IP$spc[!is.na(IP[,property]),]
Y10 = IP[,property][!is.na(IP[,property])]
IP1 = IP1[order(IP1[,property]), ]
IP1_training = IP1[c(TRUE, FALSE), ]
IP1_testing = IP1[c(FALSE, TRUE), ]
X11 = IP1_training$spc[!is.na(IP1_training[,property]),]
Y11 = IP1_training[,property][!is.na(IP1_training[,property])]
X12 = IP1_testing$spc[!is.na(IP1_testing[,property]),]
Y12 = IP1_testing[,property][!is.na(IP1_testing[,property])]

X = X10 %>% savitzkyGolay(p = 1, w = 5, m = 1)
Y = Y10
X1 = X11 %>% savitzkyGolay(p = 1, w = 5, m = 1)
Y1 = Y11
X2 = X12 %>% savitzkyGolay(p = 1, w = 5, m = 1)
Y2 = Y12
rm(X10,Y10,X11,Y11,X12,Y12)
wavelengths = as.matrix(as.numeric(colnames(X)))
matplot(x = wavelengths, y = t(X), 
        xlab = "Wavelengths (nm)",
        ylab = "Absorbance",
        type = "l", lty = 1,lwd = 2, col = "#5177A133",
        ylim = c(0, 2),
        yaxt = "n",
        cex.main = 1.2,
        cex.lab = 1.2, 
        cex.axis = 1.2,
        font.axis = 2, 
        family = "sans",
        las = 1)
y_ticks <- seq(0, 2, by = 0.5)
axis(side = 2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = 1.2)


#Second derivative calculation

#Multiplicative scatter correction

#Mean centering and scaling


####1NEAREST NEIGHBOUR####
##Choosing Dissimilarity measures##
##Training##
foo = X
foo1 = Y
#foo = X1
#foo1 = Y1
#foo = rbind(X,X1)
#rownames(foo) = as.character(seq(1, 27394, by = 1))
#foo1 = c(Y,Y1)

optimal_sel =  list(method = "opc", value = 40)
pls_tr_opc = ortho_projection(Xr = foo, Yr = foo1, method = "pls", pc_selection = optimal_sel)
pls_tr_opc
pls_comp_num = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,1]))
pls_comp_value = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,2]))
old_par = par("mfrow")
par(mfrow = c(1,1))
matplot(x = pls_comp_num, y = pls_comp_value, 
        xlab = "No of PLS components",
        ylab = "RMSD of Yr",
        type = "h", lty = 1, col = "#FF1A00CC")
title("method=pls.opc")
#Finding which k has the minimum rmse
no_components = which.min(as.numeric(pls_tr_opc$opc_evaluation[,2]))

# #PLS projection for use on new datasets
# #the pls projection matrix
# pls_tr_opc$projection_mat
# dim(pls_tr_opc$projection_mat)
# pls_projected = predict(pls_tr_opc, newdata = X2)
# dim(pls_projected)

#Comparing all the dissimilarity measures#
#PC dissimilarity with default settings (variance-based no. of components)
pcad = dissimilarity(foo, diss_method = "pca", scale = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd = dissimilarity(foo, diss_method = "pls", Yr = foo1,scale = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel = list("opc", 40)
o_pcad = dissimilarity(foo,diss_method = "pca",Yr = foo1,pc_selection = opc_sel, scale = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE)
#Correlation dissimilarity 
cd = dissimilarity(foo, diss_method = "cor", scale = TRUE)
#Moving window correlation dissimilarity 
mcd = dissimilarity(foo, diss_method = "cor", ws = 51, scale = TRUE)
#Euclidean dissimilarity 
ed = dissimilarity(foo, diss_method = "euclid", scale = TRUE)
#Cosine dissimilarity 
cosd = dissimilarity(foo, diss_method = "cosine", scale = TRUE)
#Spectral information divergence/dissimilarity 
sinfd = dissimilarity(foo, diss_method = "sid", scale = TRUE)

#Evaluations
Y_matrix = as.matrix(foo1)
ev = NULL
ev[["pcad"]] = sim_eval(pcad$dissimilarity, side_info = Y_matrix)
ev[["plsd"]] = sim_eval(plsd$dissimilarity, side_info = Y_matrix)
ev[["o_pcad"]] = sim_eval(o_pcad$dissimilarity, side_info = Y_matrix)
ev[["o_plsd"]] = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)
ev[["cd"]] = sim_eval(cd$dissimilarity, side_info = Y_matrix)
ev[["mcd"]] = sim_eval(mcd$dissimilarity, side_info = Y_matrix)
ev[["ed"]] = sim_eval(ed$dissimilarity, side_info = Y_matrix)
ev[["cosd"]] = sim_eval(cosd$dissimilarity, side_info = Y_matrix)
ev[["sinfd"]] = sim_eval(sinfd$dissimilarity, side_info = Y_matrix)

#Tabulating and plotting
statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_calib.csv", row.names = FALSE)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
par(mfrow = c(3, 3))
p = sapply(names(ev), 
            FUN = function(x, label, labs = c("Clay, %", "Clay (1-NN), %")) {
              xy = x[[label]]$first_nn[,1:2]
              plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
              title(label)
              grid()
              abline(0, 1)
              text(
                (max(xy[,1])-10), (max(xy[,2])-10),
                labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),
                pos = 1,
                col = "black", cex = 1
              )
            },
            x = ev)
par(mfrow = c(1, 1))

##Testing##
#PC dissimilarity with default settings (variance-based no. of components)
pcad = dissimilarity(Xr=foo, Xu=X2, diss_method = "pca", scale = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls", scale = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel = list("opc", 40)
o_pcad = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pca",pc_selection = opc_sel, scale = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE)
#Correlation dissimilarity 
cd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", scale = TRUE)
#Moving window correlation dissimilarity 
mcd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", ws = 51, scale = TRUE)
#Euclidean dissimilarity 
ed = dissimilarity(Xr=foo, Xu=X2, diss_method = "euclid", scale = TRUE)
#Cosine dissimilarity 
cosd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cosine", scale = TRUE)
#Spectral information divergence/dissimilarity 
sinfd = dissimilarity(Xr=foo, Xu=X2, diss_method = "sid", scale = TRUE)

#Evaluations
ev = NULL
observed_values = Y2
ev[["pcad"]]$first_nn = cbind(observed_values, foo1[apply(pcad$dissimilarity, 2, which.min)])
ev[["plsd"]]$first_nn = cbind(observed_values, foo1[apply(plsd$dissimilarity, 2, which.min)])
ev[["o_pcad"]]$first_nn = cbind(observed_values, foo1[apply(o_pcad$dissimilarity, 2, which.min)])
ev[["o_plsd"]]$first_nn = cbind(observed_values, foo1[apply(o_plsd$dissimilarity, 2, which.min)])
ev[["cd"]]$first_nn = cbind(observed_values, foo1[apply(cd$dissimilarity, 2, which.min)])
ev[["mcd"]]$first_nn = cbind(observed_values, foo1[apply(mcd$dissimilarity, 2, which.min)])
ev[["ed"]]$first_nn = cbind(observed_values, foo1[apply(ed$dissimilarity, 2, which.min)])
ev[["cosd"]]$first_nn = cbind(observed_values, foo1[apply(cosd$dissimilarity, 2, which.min)])
ev[["sinfd"]]$first_nn = cbind(observed_values, foo1[apply(sinfd$dissimilarity, 2, which.min)])

statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_valid.csv", row.names = FALSE)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
par(mfrow = c(3, 3))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("Clay, %", "Clay (1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
             title(label)
             grid()
             abline(0, 1)
             text(
               (max(xy[,1])-10), (max(xy[,2])-10),
               labels = paste("RMSD:", round(statistics_result[[label]]$RMSE,3), "\nR2:", round(statistics_result[[label]]$R_squared,3)),
               pos = 1,
               col = "black", cex = 1
             )
           },
           x = ev)
par(mfrow = c(1, 1))
rm(foo,foo1)


####2PLSR####
#foo = X
#foo1 = Y
foo = X1
foo1 = Y1
#foo = rbind(X,X1)
#rownames(foo) = as.character(seq(1, 27394, by = 1))
#foo1 = c(Y,Y1)
##Training##
pls.options(parallel = makeCluster(20, type = "PSOCK"))
#model = plsr(foo1 ~ foo, ncomp = 40, validation = "LOO")
model = plsr(foo1 ~ foo, ncomp = 40, validation = "CV", segments = 100)
stopCluster(pls.options()$parallel)
summary(model)
explvar(model)
plot(RMSEP(model), legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")
#choose the number of components with min RMSEP. 
no_components = which.min(comp_value)-1
#Since number of components are generally too large; Here we choose the number of comp with first local minimum
no_components = 18

#Plots
colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
plot(foo1, model$fitted.values[,1,no_components], xlab = "Clay %", ylab = "Clay (Global+Local PLSR) %", col = colours1[1], xlim = c(0,100), ylim = c(0,100))
grid()
abline(0, 1)
text(90, 20, labels = paste("RMSD:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$R_squared,3)),
  pos = 1,col = "black", cex = 0.8)

plot(model, plottype = "scores", comps = 1:3)
plot(model, plottype = "loadings", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Loading value")
abline(h = 0)
plot(model, plottype = "coef", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Regression coefficients")
abline(h = 0)

#Statistics calculation
statistics_result = NULL
observed_values = foo1
for (i in seq(1, 40, by = 1)) {
  predicted_values = model$fitted.values[,1,i]
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib)
write.csv(r_model_calib, "r_model_calib.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
for (i in seq(1, 40, by = 1)) {
  predicted_values = Y2_hat[,1,i] 
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
write.csv(r_model_valid, "r_model_valid.csv", row.names = FALSE)

par(mfrow = c(2,2))
matplot(r_model_calib$no_components,r_model_calib$RMSE, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_calib$no_components,r_model_calib$R_squared, xlab = "No of PLS components", ylab = "R2", type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$RMSE, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$R_squared, xlab = "No of PLS components", ylab = "R2", type = "h", lty = 1, col = "#FF1A00CC")
par(mfrow = c(1,1))
plot(Y2, Y2_hat[,1,no_components], xlab = "Clay %", ylab = "Clay (Global+Local PLSR) %", col = colours1[1], xlim = c(0,100), ylim = c(0,100))
grid()
abline(0, 1)
text(90, 20, labels = paste("N_comp:", no_components, "\nRMSD:", round(calculate_statistics(Y2,Y2_hat[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(Y2,Y2_hat[,1,no_components])$R_squared,3)),
     pos = 1,col = "black", cex = 0.8)
rm(foo,foo1)

####3KNN+PLSR####
##K-Nearest Neighbour search##
#Description of search_neighbours() function:
foo = X
foo1 = Y
#foo = X1
#foo1 = Y1
#foo = rbind(X,X1)
#rownames(foo) = as.character(seq(1, 27394, by = 1))
#foo1 = c(Y,Y1)

no_neighbor = 248
no_components = 40
knn_pc = search_neighbors(Xr = foo, Xu = X2, diss_method = "pca.nipals",k = no_neighbor)
#using PC dissimilarity with optimal selection of components
knn_opc = search_neighbors(Xr = foo, Xu = X2, diss_method = "pca.nipals",Yr = foo1,k = no_neighbor,pc_selection = list("opc", no_components),scale = TRUE)
#using PLS dissimilarity with default setting
knn_pls = search_neighbors(Xr = foo, Xu = X2, diss_method = "pls",Yr = foo1,k = no_neighbor)
#using PLS dissimilarity with optimal selection of components
knn_opls = search_neighbors(Xr = foo, Xu = X2, diss_method = "pls",Yr = foo1,k = no_neighbor,pc_selection = list("opc", no_components),scale = TRUE)
#using correlation dissimilarity
knn_c = search_neighbors(Xr = foo, Xu = X2, diss_method = "cor",k = no_neighbor, scale = TRUE)
#using moving window correlation dissimilarity
knn_mwc = search_neighbors(Xr = foo, Xu = X2, diss_method = "cor",k = no_neighbor, ws = 51, scale = TRUE)
#using Euclidean dissimilarity
knn_eu = search_neighbors(Xr = foo, Xu = X2, diss_method = "euclid",k = no_neighbor, scale = TRUE)
#using Cosine dissimilarity
knn_cos = search_neighbors(Xr = foo, Xu = X2, diss_method = "cosine",k = no_neighbor, scale = TRUE)
#using spectral information divergence dissimilarity
knn_sid = search_neighbors(Xr = foo, Xu = X2, diss_method = "sid",k = no_neighbor, scale = TRUE)
#using localized PC dissimilarity with optimal selection of components
knn_local_opc = search_neighbors(Xr = foo, Xu = X2, diss_method = "pca.nipals",Yr = foo1,k = no_neighbor,pc_selection = list("opc", no_components),scale = TRUE,.local = TRUE,pre_k = 240)
#using localized PLS dissimilarity with optimal selection of components
knn_local_opls = search_neighbors(Xr = foo, Xu = X2, diss_method = "pls",Yr = foo1,k = no_neighbor,pc_selection = list("opc", no_components),scale = TRUE,.local = TRUE,pre_k = 240)
#Printing the unique neighbours needed for all the training datasets for each method
print(length(knn_pc$unique_neighbors))
print(length(knn_opc$unique_neighbors))
print(length(knn_pls$unique_neighbors))
print(length(knn_opls$unique_neighbors))
print(length(knn_c$unique_neighbors))
print(length(knn_mwc$unique_neighbors))
print(length(knn_eu$unique_neighbors))
print(length(knn_cos$unique_neighbors))
print(length(knn_sid$unique_neighbors))
print(length(knn_local_opc$unique_neighbors))
print(length(knn_local_opls$unique_neighbors))

###Regression###
#define the regression method to be used at each neighborhood 
#Provide the optimum no of components with lowest RMSD in above NNV statistic in the dissimilarity matrix above (taken as min of RMSD)
no_components
#define the dissimilarity method from the list (pca, pls, cor, euclid, cosine, sid)
#Use Yu argument to validate the predictions directly
#Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,
#             k = seq(50, 500, by = 50), method = local_fit_pls(pls_c = no_components),
#             diss_method = "pls", pc_selection = list(method = "opc", value = 40),
#             diss_usage = "none", control = mbl_control(validation_type = "NNv"), scale = TRUE)
#Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,
#             method = local_fit_pls(pls_c = no_components),
#             k = seq(50, 500, by = 50), diss_method = knn_pc$neighbors_diss,
#             diss_usage = "none", control = mbl_control(validation_type = "NNv"), scale = TRUE)
#Y2_hat
#Plot predictions in testing data
#plot(Y2_hat, main = "Global2Local")

# #Get predictions of NNv for each combination of k and number of components
# Y2_hat$validation_results$nearest_neighbor_validation
# #collect predictions and get the indices of the best results according to nearest neighbor validation statistics
# c_val_name = "validation_results"
# c_nn_val_name = "nearest_neighbor_validation"
# c_nn_val_name = "Yu_prediction_statistics"
# Y2_hat$validation_results$Yu_prediction_statistics
# index_min = which.min(Y2_hat[[c_val_name]][[c_nn_val_name]]$rmse)
# Y2_hat_preds = get_predictions(Y2_hat)[, ..index_min]
# Y2_hat_preds = as.matrix(Y2_hat_preds)[,1]
# # Calculate statistics
# statistics_result = NULL
# statistics_result[[no_components]] <- calculate_statistics(Y2, Y2_hat_preds)


statistics_result = NULL
statistics_result1 = NULL
statistics_result2 = NULL
statistics_result3 = NULL
k1 = seq(40, 240, by = 20)
for (i in seq(1, 40, by = 1)) {
  no_components = i
  n_cores = detectCores()/2
  clust = makeCluster(n_cores)
  registerDoParallel(clust)
  Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,
               method = local_fit_pls(pls_c = no_components),
               k = seq(40, 240, by = 20), diss_method = "pca", pc_selection = list(method = "opc", value = 40),
               diss_usage = "none", control = mbl_control(validation_type = "NNv"), scale = TRUE)
  stopCluster(clust)
  #plot(Y2_hat, main = no_components)
  #collect predictions and get the indices of the best results according to nearest neighbor validation statistics
  c_val_name = "validation_results"
  c_nn_val_name = "Yu_prediction_statistics"
  #Y2_hat$validation_results$Yu_prediction_statistics
  index_min = which.min(Y2_hat[[c_val_name]][[c_nn_val_name]]$rmse)
  Y2_hat_preds = get_predictions(Y2_hat)[, ..index_min]
  Y2_hat_preds = as.matrix(Y2_hat_preds)[,1]
  print(Y2_hat$gh$projection$n_components)
  # Calculate statistics
  statistics_result[[no_components]] <- calculate_statistics(Y2, Y2_hat_preds)
  statistics_result[[no_components]]$k = k1[index_min]
  for (j in seq(1, 11, by = 1)) {
    Y2_hat_preds = get_predictions(Y2_hat)[, ..j]
    Y2_hat_preds = as.matrix(Y2_hat_preds)[,1]
    # Calculate statistics
    #statistics_result1[[as.character(no_components)]][[as.character(j)]] <- calculate_statistics(Y2, Y2_hat_preds)
    statistics_result1[[as.character(j)]] <- calculate_statistics(Y2, Y2_hat_preds)
  }
  statistics_result2 = bind_rows(statistics_result1, .id = "no_K")
  statistics_result2 = as.matrix(statistics_result2)
  statistics_result3 = rbind(statistics_result3, statistics_result2)
}

r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
write.csv(r_model_valid, "r_model_local_valid_knn_opca.csv", row.names = FALSE)
r_model_valid1 =as.matrix(statistics_result3)
print(r_model_valid1)
write.csv(r_model_valid1, "r_model_local_valid_knn_opca_all.csv", row.names = FALSE)


####Subset####
##K-Nearest Neighbour search##
foo = X
foo1 = Y
#foo = X1
#foo1 = Y1
#foo = rbind(X,X1)
#rownames(foo) = as.character(seq(1, 27394, by = 1))
#foo1 = c(Y,Y1)

no_components = 40
#PC dissimilarity with optimal selection of components
o_pcad = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pca",pc_selection = list("opc", no_components), scale = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = list("opc", no_components), scale = TRUE)

#Getting indices based on min dissimilarity values
no_neighbor = 20
unique_indices = get_lowest_indices(o_pcad$dissimilarity, no_neighbor)

##Training##
foo2 = foo[unique_indices,]
foo3 = foo1[unique_indices]
#pls.options(parallel = makeCluster(20, type = "PSOCK"))
model = plsr(foo3 ~ foo2, ncomp = 40, validation = "LOO")
explvar(model)
plot(RMSEP(model), legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")
#choose the number of components with min RMSEP. 
no_components = as.numeric(which.min(comp_value)-1)
#Since number of components are generally too large; Here we choose the number of comp with first local minimum
no_components = 18
statistics_result = NULL
foo3_hat = predict(model, newdata = foo2)
statistics_result = calculate_statistics_result(foo3, foo3_hat)
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib[no_components,])
write.csv(r_model_calib, "r_model_global_subset_20_calib.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
statistics_result = calculate_statistics_result(Y2, Y2_hat)
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
print(r_model_valid[no_components,])
write.csv(r_model_valid, "r_model_global_subset_20_valid.csv", row.names = FALSE)

#ROUGH#
print(c(min(Y), mean(Y), max(Y)))
print(c(min(Y2), mean(Y2), max(Y2)))

####PIC Band Selection####
#Data Import
IP = OSSL
rm(OSSL)
Y = as.data.frame(IP$spc[!is.na(IP$Clay),])
x = IP$Clay[!is.na(IP$Clay)]
lband = seq(1, dim(Y)[2], by = round(dim(Y)[2]/20))
lband1 = lband
lband1 = lband1[-1]
lband1 = lband1 - 1
lband1[length(lband1) + 1] <- dim(Y)[2]
#no_cores = detectCores()
no_cores = 22
#Setup cluster
clust = makeCluster(no_cores-2) #This line will take time
registerDoParallel(clust)
clusterExport(clust, c("Y","x","lband","lband1"))
all_imp_bands = foreach(i=1:20, .combine = c, .packages="NPRED") %dopar%{
  x1 = lband[i]
  x2 = lband1[i]
  py<-Y[,x1:x2]  # possible predictors
  result = stepwise.PIC(x,py)
  if (x2 == lband1[1]) {
    imp_bands = result$cpy
  } else {
    imp_bands = (x1 - 1) + result$cpy
  }}
names <- colnames(Y)
names1 = names[all_imp_bands]
names1

stopCluster(clust)
#stopImplicitCluster()
#Print the result and timing information
print(names1)





