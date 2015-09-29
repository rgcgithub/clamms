suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scales)
})

df <- read.table(
    commandArgs(trailingOnly=T)[1],
    col.names=c("sample", "cnv", "cnv.type", "copy.number",
                "exon", "in.cnv.1", "in.cnv.2", "cov",
                "max.cn.considered", "mu.dip", "sigma.dip"),
    colClasses=c("factor", "factor", "factor", "integer",
                 "factor", "factor", "factor", "numeric",
                 "integer", "numeric", "numeric")
)

df$max.cn.considered <- ifelse(df$max.cn.considered == 0, 3, -1)

sample <- levels(df$sample)[1]

for (this.cnv in levels(df$cnv)) {
    cat("Plotting ", this.cnv, "\n", sep="")
    this.cnv.df      <- filter(df, cnv == this.cnv, max.cn.considered > 0)
    this.copy.number <- this.cnv.df$copy.number[1]
    plot.max.rel.cov <- max(0.65, max((this.cnv.df$cov - this.cnv.df$mu.dip)/this.cnv.df$mu.dip), 0.5 + (this.copy.number-2)/2)

    if (nrow(this.cnv.df) > 1) {
        plot <- ggplot(this.cnv.df, aes(x = exon, y = (cov - mu.dip) / mu.dip, group=1))
        plot <- plot + scale_y_continuous(limits=c(-0.89, plot.max.rel.cov), labels=percent)
        plot <- plot + geom_ribbon(aes(ymin = -2*sigma.dip, ymax=2*sigma.dip),
                                   fill="lightgrey")
        plot <- plot + geom_ribbon(aes(ymin = -sigma.dip, ymax=sigma.dip),
                                   fill="grey")

        plot <- plot + geom_hline(yintercept=-0.5, size=1.5, colour="grey40")
        plot <- plot + geom_hline(yintercept=0.5,  size=1.5, colour="grey40")
        plot <- plot + geom_line(aes(colour = in.cnv.1), size=2)
        plot <- plot + geom_point(aes(colour = in.cnv.2), size=2)
        plot <- plot + scale_colour_manual(values=c("black", "firebrick"), guide=F)
        plot <- plot + xlab("Each tick is an exon. Grey ribbons show +/- 2\u03C3 for the diploid coverage distribution at each exon.")
        plot <- plot + ylab("Coverage relative to diploid mean")
        plot <- plot + ggtitle(paste("Sample", sample, "predicted\nto have copy number", this.copy.number, "at region", this.cnv))
        plot <- plot + theme_bw() + theme(
                    plot.title=element_text(size=8),
                    axis.title.x=element_text(size=8),
                    axis.title.y=element_text(size=8),
                    axis.text.x=element_blank(),
                    axis.text.y=element_text(size=8*5/6))

        ggsave(paste("clamms_cnv_plots/", sample, "/",
                     gsub(":", "_", gsub("-", "_", this.cnv)),
                     ".png", sep=""),
               plot=plot, units="in", width=4*1.618, height=4)
    } else {
        this.exon <- this.cnv.df$exon[1]
        cov       <- this.cnv.df$cov[1]
        mu.dip    <- this.cnv.df$mu.dip[1]
        sigma.dip <- this.cnv.df$sigma.dip[1]
        cov       <- (cov - mu.dip) / mu.dip        

        plot <- ggplot(data.frame(x=seq(-1, plot.max.rel.cov, 0.01)), aes(x=x))
        plot <- plot + scale_x_continuous(
                    limits=c(-1, plot.max.rel.cov),
                    breaks=seq(-1, plot.max.rel.cov, 0.5))
        for (k in seq(1, max(3, this.copy.number))) {
            plot <- plot + stat_function(
                        fun=dnorm,
                        arg=list(mean=0.5*(k-2), sd=sigma.dip*sqrt(k/2)))
        }
        plot <- plot + geom_point(x=cov, y=dnorm(cov, mean=0.5*(this.copy.number-2), sd=sigma.dip*sqrt(this.copy.number/2)), size=2, colour="red")
        plot <- plot + xlab("Coverage relative to diploid mean")
        plot <- plot + ylab("Probability density for mixture component")
        plot <- plot + ggtitle(paste("Sample", sample, "predicted\nto have copy number", this.copy.number, "at region", this.cnv))
        plot <- plot + theme_bw() + theme(
                    plot.title=element_text(size=6),
                    axis.title.x=element_text(size=6),
                    axis.title.y=element_text(size=6),
                    axis.text.x=element_text(size=5),
                    axis.text.y=element_text(size=5))

        ggsave(paste("clamms_cnv_plots/", sample, "/",
                     gsub(":", "_", gsub("-", "_", this.cnv)),
                     ".png", sep=""),
               plot=plot, units="in", width=2*1.618, height=2)
    }
}

