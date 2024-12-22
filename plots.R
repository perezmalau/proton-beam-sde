library("viridis")
n = 1000000

data = read.table(paste0("./test_x.txt"))
data[,2] = data[,2] / n
pdf(paste0("dose_x.pdf"))
plot(data[,1], data[,2], type="l", xlab = "Depth (cm)", ylab = "Dose", cex.lab=1.3, cex.axis=1.3, xlim=c(0,8))
dev.off()

scale_resolution = 3
data = read.table(paste0("./test_xy.txt"))
data[,3] = data[,3] / n
pdf(paste0("dose_xy.pdf"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth (cm)", ylab = "Height (cm)", cex.lab=1.3, cex.axis=1.3, cex=1, pch=15, xlim=c(0,8), ylim=c(-5,5))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol=3)
dev.off()

data = read.table(paste0("./test_z=0.txt"))
data[,3] = data[,3] / n
pdf(paste0("dose_z=0.pdf"))
offset = min(log10(data[,3]) * scale_resolution) - 1
colvec = round(log10(data[,3]) * scale_resolution - offset)
plot(data[,1], data[,2], col=viridis(max(colvec))[colvec], xlab = "Depth (cm)", ylab = "Height (cm)", cex.lab=1.3, cex.axis=1.3, cex = 1.1, pch=15, xlim=c(0,10), ylim=c(-6, 6))
legend("topleft", legend=round((min(colvec):max(colvec) + offset) / scale_resolution, 1), pch=15, col=viridis(max(colvec)), cex=1.2, ncol=3)
dev.off()
