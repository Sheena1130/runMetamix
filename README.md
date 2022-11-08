# runMetamix

The goal of this script is to implement metaMix. metaMix is a package for resolving complex metagenomic mixtures by analyzing deep sequencing data using a mixture model-based approach. 

As for more details, please check on the [metaMix package website](https://cran.r-project.org/web/packages/metaMix/). Also, you can look into [reference manual](https://cran.r-project.org/web/packages/metaMix/metaMix.pdf), [user guide](https://cran.r-project.org/web/packages/metaMix/vignettes/guide.pdf) and [related paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4565032/).


## Prerequest
Our environment is Ubuntu 20.04 (GNU/Linux 5.15.0-41-generic x86_64). Some prerequest for this implementation are as follows.

``` shell
# Install Ubuntu packages
sudo apt install openmpi-bin  # openMPI to run metaMix
sudo apt install wget  # GNU Wget to download BLAST data
sudo apt install samtools  # Samtools to get read weights
sudo apt install r-base  # R

# Install R packages
Rscript -e 'install.packages("metaMix")'
Rscript -e 'install.packages("argparse")'
```

## Preprocess
Since there are several files required for running metaMix, preprocessing is inevitable. The files prepared for running metaMix including `blast_out.txt`, `readweights.txt` and `names.dmp`, which can be obtained in the following steps. For initial users, you can use sample files offered by metaMix package simply by using the code below, and jumped to the section ***Usage of metaMix***. For other users who want to use your own data, you can go through ***Preprocess*** section to get files to run metaMix.

```shell
# Get sample files
cp /usr/local/lib/R/site-library/metaMix/extdata/blastOut_custom.tab .  # get blastOut_custom.tab
cp /usr/local/lib/R/site-library/metaMix/extdata/read_weights.tab .  # get read_weights.tab
cp /usr/local/lib/R/site-library/metaMix/extdata/names_example.dmp .  # get names_example.dmp

mv blastOut_custom.tab blast_out.txt  # rename blastOut_custom.tab as blast_out.txt
mv read_weights.tab readweights.txt  # rename read_weights.tab as readweights.txt
mv names_example.dmp names.dmp  # rename names_example.tab as names.txt
```

### Step 0. Prepare sample files
Sample files should be prepared by your own, including  `fa` and `bam`. In this implementation, we named the sample files as `sample.fa` and `sample.bam`, while `sample.fa` is contig along with reads and `sample.bam` is reads mapped to the contig.

### Step 1. Get BLAST output
The first step is to get output from BLAST, while `sample.fa` is the input file. The goal in this step is to compare a nucleotide query sequence against a protein sequence database. To get the results, you can either use [blastx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) on the [BLAST website](https://blast.ncbi.nlm.nih.gov/Blast.cgi) or use our code shown below. **Just be aware that for the BLAST site, sequences that are too long are prohibited.**

The database we use in this step is **non-redundant protein sequences (nr) from NCBI**. You can use other databases that suit your needs. Check [BLAST database](https://ftp.ncbi.nlm.nih.gov/blast/db/) to see all the options.

By using our code, you can get `blast_out.txt` as the result.

```shell
# Get BLAST bin
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar axvf ncbi-blast-2.13.0+-x64-linux.tar

# Get BLAST NCBI nr database
# [Note] Modify this part if you want to use other BLAST databases.
# [Note] This step may take a long time (can be up to 4 hours). Please be patient and stay connected to the Internet.
mkdir blast_db
cd blast_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
ls -1 nr.*.tar.gz | parallel tar axf
cd ..

# Run blastx and get BLAST output
# [Note] query should be the direction of the input fa file.
# [Note] This step may also take a long time if the sequence is long.
ncbi-blast-2.13.0+/bin/blastx -db blast_db/nr \
    -query sample.fa \
    -outfmt "6 qacc qlen sacc slen mismatch bitscore length pident evalue staxids" \
    -out blast_out.txt \
    -num_threads $(nproc)
```
If your goal is to compare a protein query sequence against a protein sequence database, you should use [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch) instead of blastx and use nucleotide collection (nr/nt) as database.

### Step 2. Get read weights
The next step of preprocessing is getting read weights from `sample.bam`. Here we use [samtools coverage](http://www.htslib.org/doc/samtools-coverage.html) as our tool and get `readweights.txt` as the result.

```shell
# Run samtools coverage to get read weights
# [Note] -H should be followed by the direction of the input bam file.
samtools coverage -H sample.bam | cut -f 1,4 > readweights.txt
```

### Step 3. Get taxonomy names file
The final step of preprocessing is to get `names.dmp`, which can simply be downloaded and extracted from NCBI.

```shell
# Download taxdump.tar.gz and get names.dmp by unzipping it
# [Note] After unzipping, not only names.dmp but also other several files will be extracted from taxdump.tar.gz. 
#        Just ignore them since they have nothing to do with the results.
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar axvf taxdump.tar.gz
```

## Usage of metaMix
After all the preprocess work, you can start to run metaMix by using `runMetamix.R`. Three files including `blast_out.txt`, `readweights.txt` and `names.dmp` will be the input of this R script. In the end, you can get `output.txt` as the final result file!

``` shell
# Run metaMix and get the final output
# [Note] Running this script may take some time for the analysis.
# [Note] If blast_out.txt comes from blastn, use 'nucl' as type instead of 'prot'.
# [Note] Threads can't be large than number of physical cores minus one (i.e. num_physical_cores - 1).
Rscript runMetamix.R \
  --reads_blast blast_out.txt \
  --reads_weights readweights.txt \
  --taxonomy_names taxdump/names.dmp \
  --output output.txt \
  --type 'prot' \
  --read_support 10 \
  --poster_prob_thr 0.9 \
  --threads 1 \
```