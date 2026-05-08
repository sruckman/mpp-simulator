# MPP Simulator

An interactive simulator for exploring multiparent population (MPP) design and QTL mapping using bulk segregant analysis (BSA).

Built for teaching how crossing design, pool size, sequencing coverage, and genetic architecture affect your ability to detect QTL.

---

## Run it in your browser (no R needed)

👉 **https://sruckman.shinyapps.io/mpp_simulator/**

---

## Run it locally in R

If you have R installed, you can run it directly from GitHub:

```r
# Install shiny if needed
install.packages("shiny")

# Run the app
shiny::runGitHub("mpp-simulator", "sruckman")
```

---

## What it simulates

- **Founders** — inbred lines whose haplotypes are mixed into the population
- **Crossing design** — fully intercrossed (all founders) or hub-and-spoke (one central founder)
- **Recombination** — generations of crossover accumulation; more generations = smaller haplotype blocks
- **Sequencing** — coverage and haplotype inference accuracy add realistic noise
- **QTL architecture** — null, 1 major QTL, 2–3 moderate QTLs, or polygenic
- **BSA pools** — cases (top N by phenotype) vs. controls (random N from unselected individuals)
- **QTL scan** — chi-square test on 2N diploid allele counts at each of 101 positions

## What the plots show

| Plot | What to look for |
|------|-----------------|
| Founder haplotypes | The pure building blocks before any recombination |
| RIL mosaic | How founders are mixed after recombination |
| Cases vs. controls | At a QTL peak, one founder color dominates in cases; controls look random |
| QTL scan | LOD score across the chromosome; red dashed line = Bonferroni threshold (LOD 3.3) |
| Founder frequency difference | Which founders are enriched in cases vs. controls at each position |
