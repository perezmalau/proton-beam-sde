library("viridis")
n = 1000000

data = read.table(paste0("./test_x.txt"))
data[,2] = data[,2] / n

pdf(paste0("bragg_peak.pdf"))
plot(data[,1], data[,2], type="l", xlab = "Depth (cm)", ylab = "Dose", cex.lab=1.3, cex.axis=1.3)
dev.off()

data = read.table(paste0("./test_xy.txt"))
data[,3] = data[,3] / n

dense_data = matrix(rep(0, round(10 * max(data[,1]) + 1) * round(10 * (max(data[,2]) - min(data[,2])) + 1)), nrow=round(10 * max(data[,1]) + 1) )
for (i in 1:dim(data)[1]) {
    dense_data[round(10 * data[i, 1]) + 1, round(10 * (data[i, 2] - min(data[,2]))) + 1] = data[i, 3]
}
col_len = 12

log_breaks = seq(log10(min(data[,3])), log10(max(dense_data)), length.out=col_len + 1)
legend_log_breaks = round(seq(log10(min(data[,3])), log10(max(dense_data)), length.out=col_len), 2)
breaks = seq(min(data[,3]), max(dense_data), length.out=col_len + 1)
legend_breaks = round(seq(min(data[,3]), max(dense_data), length.out=col_len), 2)

pdf(paste0("2d_projection_sde.pdf"))
image(x = seq(0.1, max(data[,1]) + 0.1, by=0.1), y = seq(min(data[,2]), max(data[,2]) + 0.1, by=0.1), log10(dense_data), col=viridis(col_len), breaks=log_breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="SDE")
legend("topleft", legend=legend_log_breaks, col=viridis(col_len), pch=15, bg="white", ncol=2)
dev.off()

pdf(paste0("2d_projection_sde_linear.pdf"))
image(x = seq(0.1, max(data[,1]) + 0.1, by=0.1), y = seq(min(data[,2]), max(data[,2]) + 0.1, by=0.1), dense_data, col=viridis(col_len), breaks=breaks, xlab = "Depth (cm)", ylab = "Height (cm)", main="SDE")
legend("topleft", legend=legend_breaks, col=viridis(col_len), pch=15, bg="white", ncol=2)
dev.off()