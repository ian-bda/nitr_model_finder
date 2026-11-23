# nitr_model_finder
# Find NITRs in a reference genome

Search any teleost reference genome for novel immune-type receptors (NITRs). This is the search used for rainbow trout, Atlantic salmon, and northern pike. It reports gene models, exon coordinates, and whether each locus looks complete, truncated, or like a relic. 

A NITR is an immune receptor with an immunoglobulin V domain, often an I domain, a transmembrane helix, and a cytoplasmic tail. The search starts from proteins that are already known to be NITRs and looks for relatives in the genome you supply.

## Install

Python 3.9 or newer, plus these programs on `PATH`:

```bash
conda install -c bioconda miniprot mafft hmmer samtools blast
```

That provides `miniprot`, `mafft`, `hmmbuild`, `hmmsearch`, `hmmalign`, `esl-translate`, `samtools`, `makeblastdb`, and `tblastn`.

## Run

From this directory, on a machine with enough memory (see below):

```bash
python3 find_nitrs.py \
  --genome /path/to/genome.fna.gz \
  --name my_species \
  --outdir results/my_species
```

On a cluster, the included job script asks for 8 CPUs, 80 GB of memory, and 48 hours:

```bash
mkdir -p logs
sbatch find_nitrs.sbatch /path/to/genome.fna.gz my_species results/my_species
```

`--genome` may be plain FASTA or gzipped FASTA. `--name` is a short label used in locus ids. Coordinates in the output are 1-based and inclusive on that assembly.

### Useful options

| Option | What it is |
|---|---|
| `--seeds` | Full-length NITR proteins used as queries. Default: `seeds/round1_full.faa`. |
| `--seeds-v` | Extra V-domain peptides. Default: `seeds/round1_V.faa`. |
| `--seeds-i` | I-domain peptides. The I profile is built only from this file. Default: `seeds/round1_I.faa`. |
| `--panel` | Optional full-length proteins from the species you are searching, or from a close relative. |
| `--check` | Optional proteins to compare after the search. These are never used as queries. |
| `--threads` | CPUs for miniprot and HMMER. Default: 8, or `SLURM_CPUS_PER_TASK` when set. |
| `--rebuild-rounds` | How many times to rebuild the V profile. Default: 1. |


## Steps

1. **MAFFT and HMMER `hmmbuild`.** The V peptides are aligned and turned into `NITR_V.hmm`. The I peptides are aligned separately and turned into `NITR_I.hmm`. A profile HMM is a position-by-position description of a domain, so a later search can find a similar domain after the sequence has drifted. 

2. **`miniprot`.** Each seed protein is placed on the genome, introns included. This writes the first gene models: exon coordinates, strand, and the spliced protein. Settings are `--outs 0.25`, `--outc 0.15`, and a maximum intron of 12 kb. This pass finds genes that still resemble the seeds.

3. **`esl-translate` and `hmmsearch`.** The genome is translated in all six frames. Open reading frames shorter than 60 amino acids are dropped, and anything longer than 2,000 amino acids is dropped. `hmmsearch` compares those translations to both profiles at an E-value cutoff of 1e-5. A hit is a location, not an exon model. Locations the first `miniprot` pass missed are cut out, with 4 kb of flank, and searched again with `miniprot`.

4. **The second V profile.** Proteins that passed the checks are added to the V seed set, and `hmmbuild` writes a new `NITR_V.hmm`. The checks are a GT–AG splice site, a reading frame without a stop in the middle, a V or I domain, the glycine–cysteine spacing of an Ig domain, and a methionine already in exon 1. The I profile is rebuilt from the original I peptides only, so newly found genes do not enter it. `hmmsearch` is run again on the same open reading frames. New locations go back through `miniprot`. One rebuild is the default.

5. **A second `miniprot` pass.** Complete proteins found in this genome, plus anything passed with `--panel`, are mapped onto the genome again. A distant seed is enough to notice a locus. A protein from the same fish places the leader exon more cleanly.

6. **`makeblastdb` and `tblastn`.** Published I domains, and leader peptides taken from complete models, are searched against the genome. A BLAST hit is not an exon. The script snaps each hit to the nearest GT–AG splice sites. An I hit with no gene attached to it is kept as an I-only relic.

A gene model is allowed one V and at most one I. Two V domains, or two I domains, are split into two loci. The protein is the spliced exons translated. A genomic window that still contains introns is never translated as the protein.

## Output

Everything lands in `--outdir`.

| File | Contents |
|---|---|
| `REPORT.md` | Counts of complete, truncated, relic, and review loci. |
| `tables/unique_loci.tsv` | One row per kept locus: coordinates, strand, status, domain scores. |
| `tables/exons.tsv` | One row per exon. |
| `tables/domains.tsv` | V, I, transmembrane, and tail spans on the protein. |
| `tables/unique_loci.gff3` | The same models as GFF3. |
| `proteins/unique_complete.faa` | Spliced proteins for loci called complete. Each starts with M. |
| `hmms/NITR_V.hmm` and `hmms/NITR_I.hmm` | The profiles from the last rebuild. |

`status` in the table:

- `complete` — spliced model, methionine in exon 1, a V or I domain.
- `truncated` — NITR-like, but the reading frame is broken or the start is missing.
- `relic` — a single domain on a short span, often an I domain left behind.
- `needs_review` — a non-canonical splice or some other failure of the automatic rules.

`cluster` is a label for bookkeeping. A chromosome with at least three strong V hits is called a cluster. That label is not the final cluster definition used in the trout catalog (one chromosome, gaps of at most 200 kb, at least three loci and one I domain).

Locus ids look like `my_species_mp_0001`. They are search ids. Naming a shared gene `nitr12`, a private gene `Arl_nitr1`, or a pseudogene `nitr2p` was a separate step for the trout workbook and is not part of this download.


## What this download does not do

- SignalP leader calls. The trout catalog used SignalP 6.0. Pike and salmon used SignalP 5.0. SignalP is licensed separately.
- Catalog names, presence tables, or the decision of which models to keep as pseudogenes versus drop.
- Phylogenies. Those were FAMSA alignments and IQ-TREE 2 (`-m MFP`, 1000 ultrafast bootstraps, 1000 SH-aLRT) on the V and I domains after the catalogs existed.

## Layout

```
find_nitrs.py                 command-line search
find_nitrs.sbatch             Slurm wrapper
seeds/                        published query proteins
nitrfind/gene_models.py       splice, translate, and collapse gene models
nitrfind/domain_scan.py       six-frame HMMER scan and the windows around new hits
nitrfind/domain_annotation.py leader, V, I, transmembrane, and tail labels
nitrfind/exon_repair.py       snap BLAST hits onto GT-AG exons
```

Run `find_nitrs.py`. The files in `nitrfind/` are the library it imports.
