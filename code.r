# Function to create the plot
create_circos_plot <- function() {
  library(circlize)
  
  # Read data
  chr_data <- read.csv("brassica_chromosomes.csv", stringsAsFactors = FALSE)
  homolog_pairs <- read.csv("brassica_homologs.csv", stringsAsFactors = FALSE)
  gene_density_df <- read.csv("brassica_gene_density.csv", stringsAsFactors = FALSE)
  chr_color_mapping <- read.csv("brassica_chr_colors.csv", stringsAsFactors = FALSE)
  
  chromosomes <- chr_data$chr
  chr_colors <- setNames(chr_color_mapping$color, chr_color_mapping$chr)
  genome_type <- chr_color_mapping$genome
  genome_colors <- c("A" = "#FF6B6B", "C" = "#4ECDC4")
  n_A <- sum(genome_type == "A")
  n_C <- sum(genome_type == "C")
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = c(rep(2, n_A - 1), 8, rep(2, n_C - 1), 8),
    cell.padding = c(0.02, 0, 0.02, 0),
    track.margin = c(0.005, 0.005)
  )
  
  circos.genomicInitialize(chr_data, plotType = NULL)
  
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      chr_col <- chr_colors[chr]
      circos.rect(xlim[1], 0, xlim[2], 1, col = chr_col, border = "white", lwd = 2)
      circos.text(mean(xlim), mean(ylim), chr, cex = 0.8, col = "white", font = 2,
                  facing = "bending.inside", niceFacing = TRUE)
    },
    bg.border = NA, track.height = 0.1
  )
  
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      gen_idx <- which(chr_color_mapping$chr == chr)
      gen_type <- chr_color_mapping$genome[gen_idx]
      gen_col <- genome_colors[gen_type]
      circos.rect(xlim[1], 0, xlim[2], 1, col = paste0(gen_col, "40"), border = NA)
    },
    bg.border = NA, track.height = 0.05
  )
  
  circos.genomicTrack(
    gene_density_df,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
      chr = CELL_META$sector.index
      chr_col <- chr_colors[chr]
      circos.genomicRect(region, value, col = chr_col, border = NA,
                         ytop.column = 1, ybottom = 0)
    },
    track.height = 0.1, bg.border = "#CCCCCC"
  )
  
  for(i in 1:nrow(homolog_pairs)) {
    link_col <- if(homolog_pairs$type[i] == "syntenic") "#FFA50080" else "#9370DB40"
    circos.link(homolog_pairs$chr1[i], homolog_pairs$pos1[i],
                homolog_pairs$chr2[i], homolog_pairs$pos2[i],
                col = link_col, border = NA)
  }
  
  text(-1.3, 0.5, "A genome", col = genome_colors["A"], cex = 1.2, font = 2)
  text(1.3, 0.5, "C genome", col = genome_colors["C"], cex = 1.2, font = 2)
  
  legend("topright",
         legend = c("A genome chromosomes", "C genome chromosomes", 
                    "Syntenic homologs", "Non-syntenic homologs"),
         fill = c(genome_colors["A"], genome_colors["C"], "#FFA500", "#9370DB"),
         border = c("white", "white", NA, NA), bty = "n", cex = 0.9)
  
  title("Brassica napus Genome - Homologous Genes between A and C Genomes", 
        cex.main = 1.3, font.main = 2)
  
  circos.clear()
}

# Export as JPEG (high quality)
jpeg("brassica_circos_plot.jpg", width = 3000, height = 3000, res = 300, quality = 100)
create_circos_plot()
dev.off()

# Export as PNG (lossless, high quality)
png("brassica_circos_plot.png", width = 3000, height = 3000, res = 300)
create_circos_plot()
dev.off()

# Export as PDF (vector, publication quality)
pdf("brassica_circos_plot.pdf", width = 10, height = 10)
create_circos_plot()
dev.off()

# Export as TIFF (publication quality)
tiff("brassica_circos_plot.tiff", width = 3000, height = 3000, res = 300, compression = "lzw")
create_circos_plot()
dev.off()

print("Plots exported in multiple formats!")
