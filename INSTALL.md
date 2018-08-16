# Primary processing
Primary processing -- alignment, coverage tracks, peak calling, etc -- is done in `snakemake`.

1. [Download](https://www.continuum.io) and install anaconda or miniconda, python 3.x (=not 2.x)

2. Install packages:
  ```bash
  conda install 'bwa < 0.7' samtools \
      fastx_toolkit trim-galore seqtk idr \
      ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph \
      ucsc-bigwiginfo git git-lfs htseq \
      pyBigWig weblogo twobitreader matplotlib-venn \
      wiggletools -c bioconda -c conda-forge
  ```

3. Clone the pipeline repository into `~/relmapping`:
  ```bash
  git clone https://github.com/jurgjn/relmapping.git
  ```

4. One can then submit batch jobs using `snakemake` by, at minimum, specifying `--cluster sbatch` and the number of jobs (test by adding `--dry-run`), e.g.:
  ```
  snakemake --cluster sbatch --jobs NCORES RULE
  ```

5. These two aliases (for `~/.bash_aliases`) add a few things, such as logging in sensible places :
  ```bash
  alias smj="snakemake --cluster 'sbatch -o $HOME/relmapping/logs/slurm_%j_%N.out.txt -e $HOME/relmapping/logs/slurm_%j_%N.err.txt' --use-conda --jobs "   
  alias smc="snakemake --use-conda --cores "
  ```
`smc` runs `RULE` on the head node:
  ```bash
  smc NCORES RULE
  ```
and `smj` runs `RULE` in batch mode:
  ```bash
  smj NCORES RULE
  ```
(also, add `-n` or `--dry-run` for testing a command)...

# Annotation pipeline & downstream analyses -- TODO...
