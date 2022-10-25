library("argparse")
library(metaMix)

parser <- ArgumentParser(description = 'run metamix of output file for blast')
parser$add_argument('--contig_blast', required = TRUE, help = 'path to input file')
parser$add_argument('--contig_weights', required = TRUE, help = 'path to read weight file')
parser$add_argument('--taxonomy_names', required = TRUE, help = 'path to taxonomy file')
parser$add_argument('--output', required = TRUE, help = 'path to output file')
parser$add_argument('--threads', default= 1, type = "integer", help = 'number of threads')

args <- parser$parse_args()

# read files
blastOut.file <- file.path(args$contig_blast)
read.weights <- file.path(args$contig_weights)
taxon.file <- file.path(args$taxonomy_names)

# step 1
step1 <- generative.prob(blast.output.file = blastOut.file,
                         contig.weight.file = read.weights,
                         blast.default = FALSE,
                         outDir = NULL)

# step 2
step2 <- reduce.space(step1 = step1)

# step 3
step3 <- parallel.temper(step2 = step2, noChains = args$threads)

# step 4
step4 <- bayes.model.aver(step2 = step2,
                          step3 = step3,
                          taxon.name.map = taxon.file)

# save file in txt
write.table(step4$presentSpecies.allInfo, file = args$output, quote = FALSE, sep = '\t', row.names = FALSE)
