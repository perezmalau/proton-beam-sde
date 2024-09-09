library("viridis")
n = 1000000
n_str = "1m"

data = read.table(paste0("./test_x_", n_str, ".txt"))
data[,2] = data[,2] / n
pdf(paste0("dose_x_", n_str, ".pdf"))
plot(data[,1], data[,2], type="l", xlab = "Depth", ylab = "Dose", cex.lab=1.3, cex.axis=1.3, xlim=c(0,80))
dev.off()

scale_resolution = 3
data = read.table(paste0("./test_xy_", n_str, ".txt"))
data[,3] = data[,3] / n
pdf(paste0("dose_xy_", n_str, ".pdf"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth", ylab = "Height", cex.lab=1.3, cex.axis=1.3, cex = 0.1, pch=20, xlim=c(0,80), ylim=c(-5,5))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol=3)
dev.off()
png(paste0("dose_xy_", n_str, ".png"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth", ylab = "Height", cex.lab=1.3, cex.axis=1.3, cex = 0.1, pch=20, xlim=c(0,80), ylim=c(-5,5))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol=3)
dev.off()


data = read.table(paste0("./test_z=0_", n_str, ".txt"))
data[,3] = data[,3] / n
pdf(paste0("dose_z=0_", n_str, ".pdf"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth", ylab = "Height", cex.lab=1.3, cex.axis=1.3, cex = 0.1, pch=20, ylim=c(-5, 5), xlim=c(0,80))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol = 3)
dev.off()
png(paste0("dose_z=0_", n_str, ".png"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth", ylab = "Height", cex.lab=1.3, cex.axis=1.3, cex = 0.1, pch=20, ylim=c(-5, 5), xlim=c(0,80))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol = 3)
dev.off()
