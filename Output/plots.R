library("viridis")
n = 1000000

data = read.table(paste0("./test_x.txt"))
data[,2] = data[,2] / n
geant3d = read.csv("./BP_eneDep_1E6_allProcesses.csv")

pdf(paste0("bragg_peak.pdf"))
plot(data[,1], data[,2], type="l", xlab = "Depth (cm)", ylab = "Dose", cex.lab=1.3, cex.axis=1.3, xlim=c(0,8), ylim=c(0, 5))
par(new=TRUE)
plot(geant3d[,1] / 10, geant3d[,2], type="l", xlim = c(0, 8), ylim = c(0, 5), ylab="", xlab="", xaxt="n", yaxt="n", col="red")
legend("topleft", legend=c("SDE", "Geant4"), lty=c(1,1), col=c("black", "red"), cex=1.2)
dev.off()

data = read.table(paste0("./test_xy.txt"))
data[,3] = data[,3] / n

geant2d = as.matrix(read.table("./geant2d.txt", sep=","))
ymax = round(max(data[,2]) / 0.1)
ymin = round(min(data[,2]) / 0.1)
if (ymin < -50) {
  tmp = matrix(rep(0, 100 * -round(50 + ymin)), nrow = 100)
  geant2d = cbind(tmp, geant2d)
}
if (ymax > 50) {
  tmp = matrix(rep(0, 100 * round(ymax - 50)), nrow = 100)
  geant2d = cbind(geant2d, tmp)
}

dense_data = matrix(rep(0, 100 * round(10 * (max(data[,2]) - min(data[,2])) + 1)), nrow=100)
for (i in 1:dim(data)[1]) {
    dense_data[round(10 * data[i, 1]) + 1, round(10 * (data[i, 2] - min(data[,2]))) + 1] = data[i, 3]
}
col_len = 12

log_breaks = seq(log10(min(data[,3])), log10(max(dense_data, geant2d)), length.out=col_len + 1)
legend_log_breaks = round(seq(log10(min(data[,3])), log10(max(dense_data, geant2d)), length.out=col_len), 2)
breaks = seq(min(data[,3]), max(dense_data, geant2d), length.out=col_len + 1)
legend_breaks = round(seq(min(data[,3]), max(dense_data, geant2d), length.out=col_len), 2)

pdf(paste0("2d_projection_sde.pdf"))
image(x = seq(0.1, 10, by=0.1), y = seq(min(data[,2]), max(data[,2]) + 0.1, by=0.1), log10(dense_data), col=viridis(col_len), breaks=log_breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="SDE")
legend("topright", legend=legend_log_breaks, col=viridis(col_len), pch=15, bg="white")
dev.off()

pdf(paste0("2d_projection_sde_linear.pdf"))
image(x = seq(0.1, 10, by=0.1), y = seq(min(data[,2]), max(data[,2]) + 0.1, by=0.1), dense_data, col=viridis(col_len), breaks=breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="SDE")
legend("topright", legend=legend_breaks, col=viridis(col_len), pch=15, bg="white")
dev.off()

pdf("2d_projection_geant4.pdf")
image(x = seq(0.1, 10, by=0.1), y = seq(min(data[,2]), max(data[,2]), by=0.1), log10(geant2d), col=viridis(col_len), breaks=log_breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="Geant4")
legend("topright", legend=legend_log_breaks, col=viridis(col_len), pch=15, bg="white")
dev.off()

pdf("2d_projection_geant4_linear.pdf")
image(x = seq(0.1, 10, by=0.1), y = seq(min(data[,2]), max(data[,2]), by=0.1), geant2d, col=viridis(col_len), breaks=breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="Geant4")
legend("topright", legend=legend_breaks, col=viridis(col_len), pch=15, bg="white")
dev.off()
