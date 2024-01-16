
plot.Sdescribe <- function(x, ylim = NULL, xlab = NULL, ylab = NULL, title = NULL, width = NULL, ...){



k = dim(x)[1]


resp <- trt <- NULL

df <- data.frame(trt = rownames(x), resp = x[,3])
limits <- aes(ymax = x[,5], ymin = x[,4])
out <- ggplot(df, aes(y = resp, x = trt))

if (is.null(width)) width <- 0.15 else width <- width

out <- out + geom_point() + geom_errorbar(limits, width = width, linewidth = 0.8) + theme_bw()

if (!is.null(ylim)) out <- out + ylim(ylim)

if (is.null(ylab)) out <- out + ylab("Mean and Confidence Limits") else out <- out + ylab(ylab)

if (is.null(xlab)) out <- out + xlab("Groups") else out <- out + xlab(xlab)

if (is.null(title)) out <- out + ggtitle("") else out <- out + ggtitle(title)
plot(out)


}
