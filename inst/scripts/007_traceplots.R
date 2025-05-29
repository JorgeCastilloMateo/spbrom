#################################
### Section 4 - TRACEPLOTS M5 ###
#################################

# Re-assign
model <- M5

for (i in 1:2) {
  # re-scale beta1
  model[,i]$params$beta1[,-1] <-
    sweep(model[,i]$params$beta1[,-1],  2, attr(model[,i]$x, "scaled:scale")[-1], FUN = '/')
  model[,i]$params$beta1[,1] <- model[,i]$params$beta1[,1] -
    rowSums(sweep(model[,i]$params$beta1[,-1],  2, attr(model[,i]$x, "scaled:center")[-1], FUN = '*'))
  # re-scale beta2
  model[,i]$params$beta2[,-1] <-
    sweep(model[,i]$params$beta2[,-1],  2, attr(model[,i]$x, "scaled:scale")[-1], FUN = '/')
  model[,i]$params$beta2[,1] <- model[,i]$params$beta2[,1] -
    rowSums(sweep(model[,i]$params$beta2[,-1],  2, attr(model[,i]$x, "scaled:center")[-1], FUN = '*'))
}



# CHECKING CONVERGENCE (trace-plots)
## beta1
name <- c(bquote("(Intercept) " (T[x])),
          bquote("trend " (T[x])),
          bquote("lag.tx " (T[x])),
          bquote("lag.tn " (T[x])),
          bquote("log(1+dist) " (T[x])),
          bquote("trend:lag.tx " (T[x])),
          bquote("trend:lag.tn " (T[x])),
          bquote("trend:log(1+dist) " (T[x])))
for (i in 1:8) {
  pdf(paste0("inst/img/SUPP_traceplot_",
             gsub('\\(', '', gsub('\\)', '', gsub(':', '', gsub('\\.', '',
             gsub(' ', '', colnames(model[,1]$params$beta1)[i]))))), ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$beta1[,i]), min(model[,2]$params$beta1[,i]))
  y_max <- max(max(model[,1]$params$beta1[,i]), max(model[,2]$params$beta1[,i]))
  plot(model[,1]$params$beta1[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$beta1[,i], col = "gray")
  dev.off()
}
## beta2
name <- c(bquote("(Intercept) " (T[n])),
          bquote("trend " (T[n])),
          bquote("lag.tx " (T[n])),
          bquote("lag.tn " (T[n])),
          bquote("log(1+dist) " (T[n])),
          bquote("trend:lag.tx " (T[n])),
          bquote("trend:lag.tn " (T[n])),
          bquote("trend:log(1+dist) " (T[n])))
for (i in 1:8) {
  pdf(paste0("inst/img/SUPP_traceplot_",
             gsub('\\(', '', gsub('\\)', '', gsub(':', '', gsub('\\.', '',
             gsub(' ', '', colnames(model[,1]$params$beta2)[i]))))), ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$beta2[,i]), min(model[,2]$params$beta2[,i]))
  y_max <- max(max(model[,1]$params$beta2[,i]), max(model[,2]$params$beta2[,i]))
  plot(model[,1]$params$beta2[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$beta2[,i], col = "gray")
  dev.off()
}
## a (only for M2, M4)
name <- c(expression(a["11"]), expression(a["22"]), expression(a["21"]))
for (i in 1:3) {
  pdf(paste0("inst/img/SUPP_traceplot_",
             gsub('\\(', '', gsub('\\)', '', gsub(':', '', gsub('\\.', '',
             gsub(' ', '', colnames(model[,1]$params$a)[i]))))), ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$a[,i]), min(model[,2]$params$a[,i]))
  y_max <- max(max(model[,1]$params$a[,i]), max(model[,2]$params$a[,i]))
  plot(model[,1]$params$a[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$a[,i], col = "gray")
  dev.off()
}
## a (only for M3, M5)
stations$NAME1 <-
  c("BADAJOZ", "MADRID (RETIRO)", "MALAGA", "NAVACERRADA", "SALAMANCA",
    "SAN SEBASTIAN", "TORTOSA", "VALENCIA", "ZARAGOZA", "BARCELONA (FABRA)",
    "ALBACETE", "BURGOS", "CIUDAD REAL", "CORUNA", "MURCIA", "SEVILLA",
    "SORIA", "BILBAO", "SANTIAGO", "PONFERRADA", "LEON", "LOGRONO", "ZAMORA",
    "REUS", "BARCELONA (AEROPUERTO)", "MADRID (TORREJON)", "VITORIA",
    "ALMERIA", "GIJON", "CACERES", "SANTANDER", "CASTELLON", "HUELVA",
    "LLEIDA", "MADRID (BARAJAS)", "MADRID (CUATROVIENTOS)", "MADRID (GETAFE)",
    "MORON", "VALLADOLID", "DAROCA")
for (i in 1:40) {
  # a11
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot_a11s", i, ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$a[,i]), min(model[,2]$params$a[,i]))
  y_max <- max(max(model[,1]$params$a[,i]), max(model[,2]$params$a[,i]))
  plot(model[,1]$params$a[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = bquote(a["11"](s[.(stations$NAME1[i])])))
  lines(model[,2]$params$a[,i], col = "gray")
  dev.off()
  # a22
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot_a22s", i, ".pdf"),
                  width = 6, height = 3)
  y_min <- min(min(model[,1]$params$a[,40+i]), min(model[,2]$params$a[,40+i]))
  y_max <- max(max(model[,1]$params$a[,40+i]), max(model[,2]$params$a[,40+i]))
  plot(model[,1]$params$a[,40+i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = bquote(a["11"](s[.(stations$NAME1[i])])))
  lines(model[,2]$params$a[,40+i], col = "gray")
  dev.off()
  # a21
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot_a21s", i, ".pdf"),
                  width = 6, height = 3)
  y_min <- min(min(model[,1]$params$a[,80+i]), min(model[,2]$params$a[,80+i]))
  y_max <- max(max(model[,1]$params$a[,80+i]), max(model[,2]$params$a[,80+i]))
  plot(model[,1]$params$a[,80+i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = bquote(a["11"](s[.(stations$NAME1[i])])))
  lines(model[,2]$params$a[,80+i], col = "gray")
  dev.off()
}
cat(paste0("\\includegraphics[width=.24\\textwidth]{SUPP_traceplot_a21s", 1:40, ".pdf}\n"))
## hpA (only for M3, M5)
name <- c(expression("(Intercept)" ~~~ beta[a["11"] ~ ",0"]), expression("log(1+dist)" ~~~ beta[a["11"] ~ ",1"]), expression(1 / sigma[a["11"]]^{2}), NA,
          expression("(Intercept)" ~~~ beta[a["22"] ~ ",0"]), expression("log(1+dist)" ~~~ beta[a["22"] ~ ",1"]), expression(1 / sigma[a["22"]]^{2}), NA,
          expression("(Intercept)" ~~~ beta[a["21"] ~ ",0"]), expression("log(1+dist)" ~~~ beta[a["21"] ~ ",1"]), expression(1 / sigma[a["21"]]^{2}), NA)
for (i in c(1:3,5:7,9:11)) {
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot_hp_", colnames(model[,1]$params$hpA)[i], ".pdf"),
                  width = 6, height = 3)
  y_min <- min(min(model[,1]$params$hpA[,i]), min(model[,2]$params$hpA[,i]))
  y_max <- max(max(model[,1]$params$hpA[,i]), max(model[,2]$params$hpA[,i]))
  plot(model[,1]$params$hpA[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$hpA[,i], col = "gray")
  dev.off()
}
## decay
name <- c(expression(phi), expression(phi[x]))
for (i in 1:2) {
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot_decay", i, ".pdf"),
                  width = 6, height = 3)
  y_min <- min(min(model[,1]$params$decay[,i]), min(model[,2]$params$decay[,i]))
  y_max <- max(max(model[,1]$params$decay[,i]), max(model[,2]$params$decay[,i]))
  plot(model[,1]$params$decay[,i], type = "l",
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$decay[,i], col = "gray")
  dev.off()
}
