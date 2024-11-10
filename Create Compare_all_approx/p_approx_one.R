# Загрузка пакетов
require(SPSL)
#require(gslnls)
require(minpack.lm)

# Определяем размер сетки и параметры взвешивающей функции (распределения) (b1, b2)
grid_sizes <- c(33, 65, 129, 257)
b_params <- list(c(1, 1), c(1, 2), c(2, 1))
n <- 1000  # число реализаций

# Инициализация структуры для хранения данных
results <- list()

# Цикл по значениям размера перколяционнй решетки
for (grid_size in grid_sizes) {
  
  # Создаем один элемент для записи значений
  results[[as.character(grid_size)]] <- list()
  
  # Определяем координаты сетки решетки
  xx <- seq((grid_size - 1) / 2)
  xx <- c(-rev(xx), 0, xx)
  yy <- seq(grid_size)
  
  # Внутренний цикл -- обсчет статистик при значениях параметров взвешивающего распределения (b1, b2)
  for (b in b_params) {
    b1 <- b[1]
    b2 <- b[2]
    pc <- qbeta(0.592746, b1, b2)  # оценка порога перколяции (интерактивный режим)
    cat("\nРазмер сетки:", grid_size, "b1:", b1, "b2:", b2, "pc:", pc, "\n")
    
    # Задаем граничные значени для формирования матрицы достижимости узлов решетки
    if ((b1 == 1) & (b2 == 1)) {
      pp <- c(rev(pc - 2^seq(-7, -1.00, 0.5)), pc, pc + 2^seq(-7, -1.45, 0.5))
    } else if ((b1 == 1) & (b2 == 2)) {
      pp <- c(rev(pc - 2^seq(-7, -1.01, 0.5)), pc, pc + 2^seq(-7, -1.5, 0.5))
    } else if ((b1 == 2) & (b2 == 1)) {
      pp <- c(rev(pc - 2^seq(-7, -0.5, 0.5)), pc, pc + 2^seq(-7, -2.1, 0.5))
    }
    
    # Задаем закон выбора слоев при стягивании решетки (образования кластеров)
    sst <- seq(1, grid_size) + grid_size
    trg <- grid_size^2 - 2 * grid_size + seq(1, grid_size)
    x3 <- length(trg)
    x4 <- x3 - 2
    
    # Генерация перколяционного кластера для параметров (b1, b2)
    set.seed(20241025)
    cls <- ssi20(x=grid_size, p=pc, set=sst, all=FALSE, shape=c(b1, b2))
    
    # Оцениваем Pinf для каждого массива pp
    Pinf <- rep(0, length(pp))
    cat("Вычисляем Pinf для каждого значения вероятности p в pp: ")
    for (j in seq_along(pp)) {
      cat(".")
      rfq <- fssi20(n=n, x=grid_size, p=pp[j], set=sst, all=FALSE, shape=c(b1, b2))
      Pinf[j] <- mean(rfq[trg]) * x3 / x4
    }
    cat(")\n")
    
    # Оцениваем Finf посредством cdf  в каждой точке pp[j]
    Finf <- Pinf / pbeta(pp, b1, b2)
    ppi <- seq(0, 1, 0.001)
    
    # # Оценка параметров модели по нелинейному МНК для Finf
    # fita <- gsl_nls(Finf ~ (2 / (1 + exp(-(pp - a2) / s2)) - 1) / (1 + exp(-(pp - a1) / s1)), 
    #                 start = list(a1 = pc, a2 = pc, s1 = 0.2, s2 = 0.2))
    # summary_fita <- summary(fita)
    # Fii <- predict(fita, list(pp = ppi))  
    # 
    # # Оценка параметров модели по нелинейному МНК для Pinf
    # fitb <- gsl_nls(Pinf ~ pbeta(pp, b1, b2) * (2 / (1 + exp(-(pp - a2) / s2)) - 1) / (1 + exp(-(pp - a1) / s1)), 
    #                 start = list(a1 = pc, a2 = pc, s1 = 0.2, s2 = 0.2))
    # summary_fitb <- summary(fitb)
    # Pii <- predict(fitb, list(pp = ppi))
    
    # Оценка параметров модели по нелинейному МНК для Finf
    fita <- nlsLM(Finf ~ (2 / (1 + exp(-(pp - a2) / s2)) - 1) / (1 + exp(-(pp - a1) / s1)), 
                    start = list(a1 = pc, a2 = pc, s1 = 0.03, s2 = 0.03),algorithm = "LM" )
    summary_fita <- summary(fita)
    Fii <- predict(fita, list(pp = ppi))  
    
    # Оценка параметров модели по нелинейному МНК для Pinf
    fitb <- nlsLM(Pinf ~ pbeta(pp, b1, b2) * (2 / (1 + exp(-(pp - a2) / s2)) - 1) / (1 + exp(-(pp - a1) / s1)), 
                    start = list(a1 = pc, a2 = pc, s1 = 0.03, s2 = 0.03), algorithm = "LM")
    summary_fitb <- summary(fitb)
    Pii <- predict(fitb, list(pp = ppi))
    
    
    # Сохранаяем полученные результаты
    results[[as.character(grid_size)]][[paste("b1", b1, "b2", b2, sep = "_")]] <- list(
      cls = cls,
      sst = sst,
      trg = trg,
      x3 = x3,
      x4 = x4,
      pp = pp,
      Pinf = Pinf,
      Finf = Finf,
      ppi = ppi,
      Fii = Fii,
      Pii = Pii,
      summary_fita = summary_fita,
      summary_fitb = summary_fitb
    )
  }
}

# График Pinf кривых
png("Pinf_all.png", width=14, height=10, units="in", res=1200, points=12)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(NULL, xlim=c(0.15, 1), ylim=c(0, 1), xlab="Open site probability p, 1", 
     ylab=expression("Percolation probability function "~P[infinity]*", 1"))
pch_values <- c(15, 16, 17, 18)  
cex_values <- seq(0.6, 1.2, 0.2)
ps <- c(0.05,1,0.05)

# Вывод отдельной кривой Pinf
for (i in seq_along(grid_sizes)) {
  grid_size <- grid_sizes[i]
  grid_label <- as.character(grid_size)
  
  for (j in seq_along(b_params)) {
    b <- b_params[[j]]
    b1 <- b[1]
    b2 <- b[2]
    param_label <- paste("b1=", b1, ", b2=", b2, sep="")
    
    # Извлекаем данные о Pinf кривых
    Pinf <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$Pinf
    pp <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$pp
    ppi <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$ppi
    Pii <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$Pii
    pc <- qbeta(0.592746, b1, b2)
    
    # Строим график Pinf  в зависимости от значений параметров (grid_size, b1, b2)
    points(pp, Pinf, pch=pch_values[j], cex=cex_values[i])
    #axis(1, at=seq(0.55,0.69,0.2))
    #axis(2, at=seq(0,1,0.2))
    
    lines(ppi, Pii, lwd=0.7, col="black")
    lines(pp, pbeta(pp,b1,b2), lty=4)
    abline(v=pc, col="black", lty=2)
    abline(h=c(0,1), lty=2)
  }
}

# Добавляем легенду
legend(0.12, 0.99, legend=c("S~B(1,1)", "S~B(1,2)", "S~B(2,1)"), 
       pch=pch_values, col="black", pt.cex=1, bty="n")
legend(0.12, 0.9, legend=paste("Dim(N)", grid_sizes), 
       pch=15, pt.cex=cex_values, bty="n")
dev.off()


# График Finf кривых
png("Finf_all.png", width=14, height=10, units="in", res=1200, points=12)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(NULL, xlim=c(0.15, 1), ylim=c(0, 1), xlab="Open site probability p, 1", 
     ylab=expression("Crossover function F", "1"))
pch_values <- c(15, 16, 17, 18)  
cex_values <- seq(0.6, 1.2, 0.2)

# Вывод отдельной кривой Finf
for (i in seq_along(grid_sizes)) {
  grid_size <- grid_sizes[i]
  grid_label <- as.character(grid_size)
  
  for (j in seq_along(b_params)) {
    b <- b_params[[j]]
    b1 <- b[1]
    b2 <- b[2]
    param_label <- paste("b1=", b1, ", b2=", b2, sep="")
    
    # Извлекаем данные о Finf кривых
    Finf <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$Finf
    pp <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$pp
    ppi <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$ppi
    Fii <- results[[grid_label]][[paste("b1", b1, "b2", b2, sep = "_")]]$Fii
    pc <- qbeta(0.592746, b1, b2)
    
    points(pp, Finf, pch=pch_values[j], cex=cex_values[i])
    
    lines(ppi, Fii, lwd=0.7, col="black")
    abline(v=pc, col="black", lty=2)
    abline(h=c(0,1), lty=2)
    
  }
}

# Добавляем легенду
legend(0.12, 0.99, legend=c("S~B(1,1)", "S~B(1,2)", "S~B(2,1)"), 
       pch=pch_values, col="black", pt.cex=1, bty="n")
legend(0.12, 0.9, legend=paste("Dim(N)", grid_sizes), 
       pch=15, pt.cex=cex_values, bty="n")
dev.off()


arr_1 <- sapply(results, sapply, function(l) sum(l$summary_fita$residuals^2))
arr_2 <- sapply(results, sapply, function(l) sum(l$summary_fitb$residuals^2))


