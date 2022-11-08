library("argparse")
library(metaMix)

parser <- ArgumentParser(description = 'run metamix of output file for blast')
parser$add_argument('--reads_blast', required = TRUE, help = 'path to input file')
parser$add_argument('--reads_weights', required = TRUE, help = 'path to read weight file')
parser$add_argument('--taxonomy_names', required = TRUE, help = 'path to taxonomy file')
parser$add_argument('--output', required = TRUE, help = 'path to output file')
parser$add_argument('--type', default = 'prot', required = TRUE, help = 'prot or nucl')
parser$add_argument('--read_support', default = 10, required = TRUE, help = '(for viral identification) in human clinical samples r = 10, in environmental samples r = 30')
parser$add_argument('--poster_prob_thr', default = 0.9, required = TRUE, help = 'posterior probability of presence of species threshold for reporting in the species summary')
parser$add_argument('--threads', default = 1, type = "integer", help = 'number of threads')

args <- parser$parse_args()

# read files
blastOut.file <- file.path(args$reads_blast)
read.weights <- file.path(args$reads_weights)
taxon.file <- file.path(args$taxonomy_names)

# step 1
if (args$type == 'prot') {
   step1 <- generative.prob(blast.output.file = blastOut.file,
                         contig.weight.file = read.weights,
                         blast.default = FALSE,
                         outDir = NULL)
} else {
   step1 <- generative.prob.nucl(blast.output.file = blastOut.file,
                         contig.weight.file = read.weights,
                         blast.default = FALSE,
                         outDir = NULL)
}


# step 2
step2 <- reduce.space(step1 = step1)

# step 3
if (args$type == 'prot') {
   step3 <- parallel.temper(step2 = step2, readSupport = args$read_support, noChains = args$threads)
} else {
   step3 <- parallel.temper.nucl(step2 = step2, readSupport = args$read_support, noChains = args$threads)
}

# step 4
step4 <- bayes.model.aver(step2 = step2,
                          step3 = step3,
                          taxon.name.map = taxon.file,
                          poster.prob.thr = args$poster_prob_thr)

# save file in txt
write.table(step4$presentSpecies.allInfo, file = args$output, quote = FALSE, sep = '\t', row.names = FALSE)
