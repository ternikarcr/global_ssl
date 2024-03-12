packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", 
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
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

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
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

setwd("D:/Academics/PhD/SEM 11/Paper_3_Global2Local/Working_files/Refined/")
setwd("F:/CRT/Refined/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
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
new_wavs = seq(502, 2500, by = 2)
foo1$spc = foo2 %>% 
  resample(wav = old_wavs, new.wav = new_wavs)
kar = foo1
rm(foo,foo1,foo2)

#Berambadi Data
foo = read_excel("Berambadi_IP.xlsx", sheet = "Sheet1")
foo1 = as.data.frame(foo[,1:3])
foo2 = -log10(data.matrix(foo[,4:length(foo)]))
rownames(foo2) = as.character(seq(1, 275, by = 1))
colnames(foo2) = as.character(seq(350, 2500, by = 1))
#pre-process the spectra: resample it to a resolution of 2 nm 
old_wavs = foo2 %>% colnames() %>% as.numeric()
new_wavs = seq(502, 2500, by = 2)
foo1$spc = foo2 %>% 
  resample(wav = old_wavs, new.wav = new_wavs)
Ber = foo1
rm(foo,foo1,foo2)

#Belgium Data
foo = read_csv("Belgium.csv")
foo1 = as.data.frame(foo[,701:703])
foo2 = -log10(data.matrix(foo[,1:700]))
rownames(foo2) = as.character(seq(1, 825, by = 1))
colnames(foo2) = as.character(seq(1100, 2498, by = 2))
foo1$SOC[foo1$SOC == -9999] <- NA
foo1$CEC[foo1$CEC == -9999] <- NA
foo1$N[foo1$N == -9999] <- NA
foo1$spc = foo2 
Belgium = foo1
rm(foo,foo1,foo2)

#Data provision
IP = OSSL
IP = kar
IP = Ber
IP = Belgium
#Component selection
property = "N"
summary(IP[,property])

imp_no_comp = list()
imp_no_samples = list()
IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For clay, silt, sand the values were made to selected in range of 0-100
#IP_Y1[IP_Y1 < 0 | IP_Y1 > 100] <- NA
#For Avail_P the values were made to selected in range of 0-200
#IP_Y1[IP_Y1 < 0.5 | IP_Y1 > 200] <- NA
#IP_Y = IP_Y1[!is.na(IP_Y1)]
#IP_X = IP_X1[!is.na(IP_Y1),]
#summary(IP_Y)
#IP_X_O = IP_X[order(IP_Y), ]
#IP_Y_O = IP_Y[order(IP_Y)]
#Split 10:1
# X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 5:1
# X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 2:1
# X1 = IP_X_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 1:1
# X1 = IP_X_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#For Karnataka, Berambadi, Belgium
#Split 5:1
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
#Split 2:1
X1 = IP_X_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
#Split 1:1
X1 = IP_X_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

####1NEAREST NEIGHBOUR####
##Choosing Dissimilarity measures##
##Training##
foo = X1
foo1 = Y1

#Optional
optimal_sel =  list(method = "opc", value = min(length(foo1)-4,40))
#optimal_sel =  list(method = "manual", value = 27)
#no_components = 27
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
print(paste("Number of components:", no_components))
imp_no_comp = append(imp_no_comp, no_components)

#Comparing all the dissimilarity measures#
#opc_sel = list("opc", 40)
opc_sel = list("manual", no_components)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE)
o_plsd_tr = o_plsd

#Evaluations
Y_matrix = as.matrix(foo1)
ev = NULL
ev[["o_plsd"]] = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)

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
par(mfrow = c(1, 1))
xy = ev[[label]]$first_nn[,1:2]
plot(xy[,1], xy[,2], xlab = property, ylab = paste0(property," predicted") , col = colours1[1])
grid()
abline(0, 1)
text((max(xy[,1])-10), (max(xy[,2])-10), labels = paste("RMSD:", round(as.numeric(statistics_result[[label]][2]),3), "\nR2:", round(as.numeric(statistics_result[[label]][1]),3)),
               pos = 1,col = "black", cex = 1)
par(mfrow = c(1, 1))

##Testing##
#PC dissimilarity with optimal selection of components
#opc_sel = list("opc", 40)
opc_sel = list("manual", no_components)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE)
o_plsd_test = o_plsd

#Evaluations
ev = NULL
observed_values = Y2
ev[["o_plsd"]]$first_nn = cbind(observed_values, foo1[apply(o_plsd$dissimilarity, 2, which.min)])

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
par(mfrow = c(1, 1))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("Clay", "Clay predicted")) {
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
#rm(foo,foo1)


####2PLSR####
##Training##
#pls.options(parallel = makeCluster(20, type = "PSOCK"))
if (length(foo1) < 1000) {
  model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,40), validation = "LOO")  
} else {
  model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,40), validation = "CV", segments = 1000)
}
#stopCluster(pls.options()$parallel)
summary(model)
explvar(model)
plot(RMSEP(model), legendpos = "topright")
#plot(model$validation$PRESS[1,], legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")

#choose the number of components with min RMSEP. 
no_components = which.min(comp_value)-1
#Since number of components are generally too large; Here we choose the number of comp with first local minimum
#no_components = 20
no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
#no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
#no_components = min(which.min(comp_value)-1,no_components_onesigma, no_components_permut)
no_components = min(which.min(comp_value)-1,no_components_onesigma)
imp_no_comp = append(imp_no_comp, no_components)
imp_no_samples = append(imp_no_samples, length(foo1))

#Plots
colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
plot(foo1, model$fitted.values[,1,no_components], xlab = property, ylab = paste0(property," predicted"), col = colours1[1], xlim = c(0,100), ylim = c(0,100))
grid()
abline(0, 1)
text(90, 30, labels = paste("N_comp:", no_components, "\nRMSD:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$R_squared,3)),
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
for (i in seq(1, min(length(foo1)-1,40), by = 1)) {
  predicted_values = model$fitted.values[,1,i]
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib)
print(r_model_calib[no_components,])
write.csv(r_model_calib, "r_model_calib.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
for (i in seq(1, min(length(foo1)-4,40), by = 1)) {
  predicted_values = Y2_hat[,1,i] 
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
print(r_model_valid[no_components,])
write.csv(r_model_valid, "r_model_valid.csv", row.names = FALSE)

par(mfrow = c(2,2))
matplot(r_model_calib$no_components,r_model_calib$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_calib$no_components,r_model_calib$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")
par(mfrow = c(1,1))
plot(Y2, Y2_hat[,1,no_components], xlab = property, ylab = paste0(property," predicted"), col = colours1[1], xlim = c(0,100), ylim = c(0,100))
grid()
abline(0, 1)
text(90, 30, labels = paste("N_comp:", no_components, "\nRMSD:", round(calculate_statistics(Y2,Y2_hat[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(Y2,Y2_hat[,1,no_components])$R_squared,3)),
     pos = 1,col = "black", cex = 0.8)
#rm(foo,foo1)


####Subset####
##K-Nearest Neighbour search##
#Getting indices based on min dissimilarity values
for (no_neighbor in list(1,2,5,10)) {
  unique_indices = get_lowest_indices(o_plsd$dissimilarity, no_neighbor)
  print(length(unique_indices))
  imp_no_samples = append(imp_no_samples, length(unique_indices))
  
  ##Training##
  foo2 = foo[unique_indices,]
  foo3 = foo1[unique_indices]
  #pls.options(parallel = makeCluster(20, type = "PSOCK"))
  if (length(foo3) < 1000) {
    model = plsr(foo3 ~ foo2, ncomp = min(length(unique_indices)-4,40), validation = "LOO")  
  } else {
    model = plsr(foo3 ~ foo2, ncomp = min(length(unique_indices)-4,40), validation = "CV", segments = 1000)
  }
  explvar(model)
  plot(RMSEP(model), legendpos = "topright")
  par(mfrow = c(1,1))
  comp_value = RMSEP(model)$val[1,1,]
  matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")
  
  #choose the number of components with min RMSEP. 
  no_components = which.min(comp_value)-1
  #Since number of components are generally too large; Here we choose the number of comp with first local minimum
  #no_components = 20
  no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
  #no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
  #no_components = min(which.min(comp_value)-1,no_components_onesigma, no_components_permut)
  no_components = min(which.min(comp_value)-1,no_components_onesigma)
  no_components = max(2,no_components_onesigma)
  imp_no_comp = append(imp_no_comp, no_components)
  
  statistics_result = NULL
  foo3_hat = predict(model, newdata = foo2)
  statistics_result = calculate_statistics_result(foo3, foo3_hat, min(length(unique_indices)-4,40))
  r_model_calib = bind_rows(statistics_result, .id = "no_components")
  print(r_model_calib[no_components,])
  write.csv(r_model_calib, paste0("r_model_nn_",no_neighbor,"_calib.csv"), row.names = FALSE)
  
  ##Testing##
  statistics_result = NULL
  observed_values = Y2
  Y2_hat = predict(model, newdata = X2)
  statistics_result = calculate_statistics_result(Y2, Y2_hat, min(length(unique_indices)-4,40))
  r_model_valid = bind_rows(statistics_result, .id = "no_components")
  print(r_model_valid)
  print(r_model_valid[no_components,])
  write.csv(r_model_valid, paste0("r_model_nn_",no_neighbor,"_valid.csv"), row.names = FALSE)
}


####Uniform model####
IP_X = IP$spc[!is.na(IP[,property]),]
IP_Y = IP[,property][!is.na(IP[,property])]
IP_X_O = IP_X[order(IP_Y), ]
IP_Y_O = IP_Y[order(IP_Y)]

#num_segments = 40
num_segments = 200
segments = cut(IP_Y_O, breaks = seq(0, max(IP_Y_O), length.out = num_segments+1), include.lowest = TRUE)
grouped_Y = split(IP_Y_O, segments)
#selected_Y = sapply(grouped_Y, head, n = 3)
selected_Y = sapply(grouped_Y, head, n = 1)
indices = unlist(lapply(selected_Y, function(sample) match(sample,IP_Y_O)))
length(indices)
length(unique(indices))
foo2 = IP_X_O[indices,]
foo3 = IP_Y_O[indices]
hist(IP_Y_O)
hist(foo3)
rm(IP_X,IP_Y,IP_X_O,IP_Y_O)
rm(selected_Y)

#pls.options(parallel = makeCluster(20, type = "PSOCK"))
model = plsr(foo3 ~ foo2, ncomp = 40, validation = "LOO")
explvar(model)
plot(RMSEP(model), legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Yr", type = "h", lty = 1, col = "#FF1A00CC")

#choose the number of components with min RMSEP. 
#Since number of components are generally too large; Here we choose the number of comp with first local minimum
#no_components = 20
no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
#no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
#no_components = min(which.min(comp_value)-1,no_components_onesigma, no_components_permut)
no_components = min(which.min(comp_value)-1,no_components_onesigma)
imp_no_comp = append(imp_no_comp, no_components)
imp_no_samples = append(imp_no_samples, length(foo3))

statistics_result = NULL
foo3_hat = predict(model, newdata = foo2)
statistics_result = calculate_statistics_result(foo3, foo3_hat, min(length(unique_indices)-4,40))
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib[no_components,])
write.csv(r_model_calib, "r_model_uni_calib.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
statistics_result = calculate_statistics_result(Y2, Y2_hat, min(length(unique_indices)-4,40))
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
print(r_model_valid[no_components,])
write.csv(r_model_valid, "r_model_uni_valid.csv", row.names = FALSE)


####ROUGH####

####3KNN+PLSR####
##K-Nearest Neighbour search##
statistics_result = NULL
statistics_result1 = NULL
statistics_result2 = NULL
statistics_result3 = NULL
k1 = seq(40, 240, by = 20)

no_components
Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,
             method = local_fit_pls(pls_c = no_components),
             k = seq(40, 240, by = 20), diss_method = o_plsd$dissimilarity,
             diss_usage = "none", control = mbl_control(validation_type = "NNv"), scale = TRUE)
#plot(Y2_hat, main = no_components)
#collect predictions and get the indices of the best results according to nearest neighbor validation statistics
c_val_name = "validation_results"
c_nn_val_name = "Yu_prediction_statistics"
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
  statistics_result1[[as.character(j)]] <- calculate_statistics(Y2, Y2_hat_preds)
}
statistics_result2 = bind_rows(statistics_result1, .id = "no_K")
statistics_result2 = as.matrix(statistics_result2)
statistics_result3 = rbind(statistics_result3, statistics_result2)


r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
write.csv(r_model_valid, "r_model_global_valid_knn_opls.csv", row.names = FALSE)
r_model_valid1 =as.matrix(statistics_result3)
print(r_model_valid1)
write.csv(r_model_valid1, "r_model_global_valid_knn_opls_all.csv", row.names = FALSE)

####ROUGH####
#Transformations
# Example data
hist(IP[, property][IP[, property] >= 1 & IP[, property] <= 200])
foo_IP = IP[, property][IP[, property] >= 1 & IP[, property] <= 200]
foo_IP1 = foo_IP[!is.na(foo_IP)]
summary(foo_IP1)
original_data <- foo_IP1
summary(original_data)

# Log transformation
log_transformed <- log(original_data)
# Square root transformation
sqrt_transformed <- sqrt(original_data)
# Box-Cox transformation
library(MASS)
bc_result <- boxcox(original_data ~ 1)
# Extract lambda value
lambda <- bc_result$x[which.max(bc_result$y)]
# Apply Box-Cox transformation
if (lambda == 0) {boxcox_transformed <- log(original_data)} else {
  boxcox_transformed <- (original_data^lambda - 1) / lambda
}

# Inverse transformation
inverse_transformed <- 1 / original_data
# Rank transformation
ranked_data <- rank(original_data)
rank_transformed <- qnorm((ranked_data - 0.5) / length(original_data))
# Quantile transformation
quantile_transformed <- qnorm(ppoints(original_data))
# Output
transformations <- list(
  Log = log_transformed,
  Square_Root = sqrt_transformed,
  Box_Cox = boxcox_transformed,
  Inverse = inverse_transformed,
  Rank = rank_transformed,
  Quantile = quantile_transformed
)
# Plot histograms of transformed data
par(mfrow = c(4, 2))
hist(original_data, main = "Original Data", breaks = 20)
for (i in 1:length(transformations)) {
  hist(transformations[[i]], main = names(transformations)[i], breaks = 20)
}

par(mfrow = c(1, 1))
# Plot histograms of transformed data in one figure
hist(original_data, main = "Transformed Data Comparison", breaks = 20, prob = TRUE, col = "lightgray", xlab = "Value", ylab = "Density")

# Plot the histograms of transformed data, each with a different color
for (i in 1:length(transformations)) {
  hist(transformations[[i]], prob = TRUE, add = TRUE, col = rainbow(length(transformations))[i])
}

# Quantile transformation
quantile_transformed <- qnorm(ppoints(original_data))
# Inverse quantile transformation from quantile-transformed data to original data
inverse_transformed_data <- qnorm(quantile_transformed)
summary(original_data)
summary(quantile_transformed)
summary(inverse_transformed_data)
par(mfrow = c(2, 1))
hist(quantile_transformed, main = "Transformed Data Comparison", breaks = 20, prob = TRUE, col = "lightgray", xlab = "Value", ylab = "Density")
hist(inverse_transformed_data, main = "Transformed Data Comparison", breaks = 20, prob = TRUE, col = "lightgray", xlab = "Value", ylab = "Density")
