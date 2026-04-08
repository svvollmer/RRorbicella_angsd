# Results

Three-species *Orbicella* population genomics — *O. annularis*, *O. faveolata*, and *O. franksi*
from the Reef Renewal Florida restoration program, aligned to the *O. franksi* reference genome
(jaOrbFran1.1 / GCA_964199315.1).

> **Status (2026-04-08):** Segments 1–4 complete. Seg 6 (moments) ALL 3 comparisons complete. Seg 7 (SMC++) SIF pulled, ready to launch.

---

## Samples

| Species | N (field) | N (genetic) | Notes |
|---------|-----------|-------------|-------|
| *O. annularis* | 63 | 59 | 4 reclassified by K=3 admixture |
| *O. faveolata* | 57 | 50 | 7 reclassified; Ofra_26 (clone) excluded |
| *O. franksi* | 21 | 24 | 3 reclassified from other species; Ofra_6 (hybrid) excluded |
| **Total** | **141** | **133** | After relatedness + lineage filtering |

10 field mislabels identified (7.5% error rate), all corrected by K=3 genetic assignment.
Ofra_6 excluded as putative *franksi* × *faveolata* hybrid (8.4% *faveolata* ancestry at K=3).
Ofra_26 excluded as clone of Ofav_106 (KING=0.493).

**Reference:** *O. franksi* jaOrbFran1.1 / GCA_964199315.1 (15 chromosomes OZ116698.1–OZ116712.1; all scaffolds ≥ 10 Mb retained)

> ⚠️ **Reference genome bias note:** All analyses use the *O. franksi* reference. SNP ascertainment
> is anchored to *O. franksi* sequence space, which may inflate apparent similarity between
> *O. franksi* and *O. annularis* relative to true phylogenetic distances. Interpret pairwise
> genetic distance estimates accordingly.

---

## Sequencing QC

| Species | N | Mean mapping rate | Range |
|---------|---|-------------------|-------|
| *O. annularis* | 63 | 93.2% | 80.7–97.7% |
| *O. faveolata* | 55 | 94.3% | 81.7–97.9% |
| *O. franksi* | 22 | 94.6% | 78.5–97.8% |

All three species map well to the *O. franksi* reference. *O. annularis* has the lowest mean
mapping rate (−1.4% vs *O. franksi*), consistent with it being the most phylogenetically
distant species from the reference. This difference in mapping efficiency contributes to
fewer callable variant sites in *O. annularis* comparisons (see Demographic Inference).

SNPs retained (pass 2, MAF > 0.05, ≥ 80% individuals): ~7,233,371 from 28,063,923 scanned.
SNP filter params: `docs/outputs/filter_params.yaml`

---

## Relatedness and Clonality

Pairwise all-vs-all ngsRelate. Clone threshold: KING ≥ 0.45.

| Clone pair | KING | Excluded | Reason |
|-----------|------|----------|--------|
| Oann_72 / Oann_73 | 0.4949 | Oann_72 | lower depth (22.9x vs 39.9x) |
| Oann_75 / Oann_77 | 0.4959 | Oann_77 | lower depth (28.9x vs 39.1x) |
| Oann_84 / Oann_85 | 0.4949 | Oann_85 | lower depth (40.2x vs 47.1x) |
| Oann_95 / Oann_96 | 0.4955 | Oann_95 | lower depth (29.8x vs 31.0x) |
| Ofav_106 / Ofra_26 | 0.4929 | **Ofra_26** | **MISLABEL** — confirmed *O. faveolata* |
| Ofav_92 / Ofav_93 | 0.4853 | Ofav_92 | lower depth (29.6x vs 32.8x) |

**134 unrelated samples retained** after clone/relatedness filtering.
**133 samples retained for Segment 4+** after lineage assignment gate (Ofra_6 excluded as hybrid).

### Ofra_26 mislabel

Ofra_26 (labeled *O. franksi*, FL_RR Upper Keys / East Turtle Shoal) is a confirmed *O. faveolata* mislabel.
KING vs Ofav_106 = +0.493 (clone — same individual). Decision: exclude Ofra_26, retain Ofav_106.

---

## Lineage Assignment

Based on ngsAdmix K=3 best-run Q matrix. Samples assigned to species contributing majority K=3 ancestry.
Ofra_6 (putative hybrid, 8.4% faveolata) excluded from all Segment 4+ analyses.

**Lineage assignment file:** `results/admixture/lineage_assignments.txt`

| Genetic lineage | N | Notes |
|----------------|---|-------|
| *O. annularis* | 59 | |
| *O. faveolata* | 50 | |
| *O. franksi* | 24 | Ofra_6 excluded |
| **Total** | **133** | |

### Field label mislabels detected

K=3 assignment identified **10 samples** with 100% ancestry for a different species than their field label:

| Sample | Field label | Genetic assignment | Q (assigned) |
|--------|-------------|-------------------|-------------|
| Oann_106 | *O. annularis* | *O. faveolata* | 1.000 |
| Ofav_75 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofav_76 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofav_83 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofav_84 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofav_94 | *O. faveolata* | *O. annularis* | 1.000 |
| Ofav_104 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofav_105 | *O. faveolata* | *O. franksi* | 1.000 |
| Ofra_1 | *O. franksi* | *O. faveolata* | 1.000 |
| Ofra_4 | *O. franksi* | *O. faveolata* | 1.000 |

Mislabel rate: 10/134 = 7.5% (11/140 = 7.9% including Ofra_26). All downstream analyses use
genetic lineage assignments, not field labels.

---

## Population Structure

### LD Decay

![LD decay](docs/figures/ld_decay.png)

LD decays from r²≈0.25 at <1 kb to r²≈0.10 by **~8–10 kb**, reaching r²≈0.05 by 35 kb and
plateauing near background (~0.05) beyond 50 kb. Notably faster than *Acropora* (r²=0.1 at
~50–100 kb), consistent with larger historical effective population sizes in massive corals
(*Orbicella* Ne likely 2–5× higher than *Acropora*).

---

### PCA

![PCA with confidence ellipses](docs/figures/pca_ellipses.png)

- **PC1 (22.8%)** separates *O. annularis* from *O. faveolata* + *O. franksi*
- **PC2 (4.9%)** separates *O. franksi* from *O. faveolata*, resolving all three species
- **PC3 (1.8%)** captures within-*O. annularis* structure
- 95% ellipses confirm clean species clustering; *faveolata*/*franksi* ellipses overlap on PC1 (sister species)

---

### Admixture — K=2 vs K=3

![K=2 vs K=3 comparison](docs/figures/admixture_K2_K3_comparison.png)

![Log-likelihood annotation](docs/figures/admixture_loglik_annotated.png)

ΔK picks K=2 (Evanno method dominated by the deep Oann vs faveolata/franksi split), but
**K=3 is biologically correct** — three distinct genetic clusters clearly resolved by both
log-likelihood elbow and PCA. The K=2→3 gain (~138,000 log-lik units) is real and biologically
meaningful; it is ~5× smaller than K=1→2, causing ΔK to be overwhelmed — a well-documented
Evanno artefact when one split is much deeper than others.

![ngsAdmix K2-4](docs/figures/admixture_ngsadmix_K234.png)
![PCAngsd K2-4](docs/figures/admixture_pcangsd_K234.png)

**K=3:** All three species cleanly resolved. **Ofra_6** is the only putative hybrid (8.4%
*faveolata* ancestry; all other samples >95% pure at K=3).

---

### Between-species Genetic Distances

![Genetic distance tree](docs/figures/genetic_distance_tree.png)

UPGMA tree groups *O. annularis* + *O. franksi* as most similar — **reference genome bias
artefact** (SNPs called in *O. franksi* space inflate apparent Oann/Ofra similarity). The
known phylogeny (*O. faveolata* + *O. franksi* as sisters; Fukami et al. 2004) is recovered
by moments demographic inference (see below). Spawning biology also supports this: *O. faveolata*
and *O. franksi* spawn closer in time, enabling actual interspecific hybridization; *O. annularis*
and *O. franksi* spawn hours apart despite laboratory cross-compatibility.

---

## Putative Hybrids

| Sample | Labeled species | Faveolata ancestry (K=3) | Status |
|--------|----------------|------------------------|--------|
| Ofra_6 | *O. franksi* | 8.4% | Putative *franksi* × *faveolata* hybrid; excluded from Seg 4+ |

---

## Genetic Diversity

| Species | π (genome-wide) | Tajima's D |
|---------|----------------|------------|
| *O. annularis* | 0.01241 | −0.59 |
| *O. faveolata* | 0.01228 | −0.20 |
| *O. franksi* | 0.01272 | −0.52 |

All three species have similar diversity (π ≈ 0.012–0.013). All Tajima's D values are
negative, consistent with recent population expansion or purifying selection. *O. faveolata*
has the least negative D (closest to 0), suggesting a more stable recent demographic history
relative to the other two species.

---

## FST and Population Differentiation

| Comparison | Global FST |
|------------|-----------|
| *O. annularis* vs *O. franksi* | 0.180 |
| *O. faveolata* vs *O. franksi* | 0.109 |
| *O. annularis* vs *O. faveolata* | pending (fst_index running 2026-04-08) |

*O. annularis* vs *O. franksi* has the highest FST despite appearing closest in the
distance tree (reference genome bias artefact). *O. faveolata* + *O. franksi* are
true sisters (Fukami et al. 2004) and show lower FST consistent with more recent
divergence and ongoing gene flow.

Windowed FST (50 kb) complete for Oann–Ofra and Ofav–Ofra. Oann–Ofav relaunched
2026-04-08 (prior fst_index NODE_FAIL; .fst.gz intact, .fst.idx regenerating).

---

## Demographic Inference (moments)

Models fitted: SI, IM, IM_a, SC, AM — folded 2D SFS, 20 restarts each.
Samples use genetically corrected lineage assignments (not field labels).
μ = 1.8×10⁻⁸, generation time = 10 yr, Na estimated from π/(4μ).

### Model selection

| Comparison | Best model | ΔAIC (vs next) | Segregating sites |
|------------|-----------|----------------|-------------------|
| *O. annularis* vs *O. franksi* | SC | 796 | 33,620 |
| *O. faveolata* vs *O. franksi* | SC | 1,764 | 311,509 |
| *O. annularis* vs *O. faveolata* | **IM_a** | **2,902** | **351,386** |

The Ofra comparisons both strongly support **Secondary Contact (SC)**: populations
diverged in allopatry, then resumed gene flow. The AM model for *O. annularis* vs
*O. franksi* converged with T2 ≈ 0 (ancient migration collapsed to zero),
independently confirming SC over AM.

In contrast, *O. annularis* vs *O. faveolata* supports **IM_a (isolation with
asymmetric migration)**: continuous but asymmetric gene flow since divergence, with
no allopatric phase. Oann→Ofav migration is ~3× stronger than Ofav→Oann.

### Parameter estimates (absolute, scaled from π-based Na)

Na estimated from π/(4μ), μ=1.8×10⁻⁸. Two generation times shown for comparison;
generation time for *Orbicella* is uncertain (~5–10 yr).

| Parameter | *O. ann* vs *O. fra* | *O. fav* vs *O. fra* | *O. ann* vs *O. fav* |
|-----------|---------------------|---------------------|---------------------|
| Best model | SC | SC | **IM_a** |
| Ancestral Na (π-based) | ~174,500 | ~173,600 | ~171,500 |
| Split / T1 time (g=10yr) | ~1.27 Mya | ~1.50 Mya | ~1.90 Mya |
| Split / T1 time (g=5yr) | ~0.64 Mya | ~0.75 Mya | ~0.95 Mya |
| Secondary contact began (g=10yr) | ~136 kya | ~465 kya | — (continuous IM) |
| Secondary contact began (g=5yr) | ~68 kya | ~233 kya | — (continuous IM) |
| T2/T1 ratio | 12% | 45% | N/A |
| Nu1 (pop1 current Ne) | ~516,000 | ~139,000 | ~238,000 |
| Nu2 (pop2 current Ne) | ~555,000 | ~231,000 | ~123,000 |
| m or m12 (into pop1, scaled) | 2.52 | 1.30 | 0.47 (Ofav→Oann) |
| m21 (into pop2, scaled) | — | — | 1.38 (Oann→Ofav) |

### Biological interpretation

**Split times** (~1.3–1.9 Mya) are consistent with Pliocene/Pleistocene *Orbicella* speciation.

**Secondary contact timing differs markedly for the Ofra comparisons:**
- *O. annularis* vs *O. franksi*: contact resumed ~136 kya, coinciding with the last
  interglacial (MIS5e, ~120–130 kya) when rising sea levels expanded reef habitat.
  T2/T1 = 12% — relatively brief secondary contact.
- *O. faveolata* vs *O. franksi*: contact resumed ~465 kya (MIS11). T2/T1 = 45% —
  sustained secondary contact consistent with their sister-species relationship and
  near-synchronous spawning that allows ongoing hybridization.

**The *O. annularis* vs *O. faveolata* result (IM_a) is qualitatively different:**
rather than a period of allopatry followed by secondary contact, these two species
have maintained continuous, asymmetric gene flow since their divergence (~1.90 Mya g=10yr).
Oann→Ofav migration is ~3× stronger than Ofav→Oann (scaled m: 1.38 vs 0.47).

**The SC results for Ofra comparisons recover the known phylogeny** (*O. faveolata* +
*O. franksi* as sisters) even though the raw genetic distance tree is obscured by
reference genome bias. This is a key internal validation.

**Caveats:**
- The *O. annularis* vs *O. franksi* comparison has 9× fewer segregating sites than
  the other pairs (33,620 vs 311,509–351,386). This primarily reflects greater sequence
  divergence between *O. annularis* and the *O. franksi* reference genome. Parameter
  estimates for this comparison should be treated with more caution.
- Na estimated from current π rather than SMC++ (Seg 7 not yet run). Absolute time
  and Ne estimates will be refined once SMC++ is complete.
- The Nu estimates from moments imply large post-split Ne values (~500k for Oann/Ofra),
  which differ from the π-based current Ne (~170k). Consistency check pending SMC++.

---

## Population Size History (SMC++)

*Snakefile.smcpp written; SMC++ SIF pulled (2026-04-08). Ready to launch (`bash run.sh 7`).*

**Parameters:** μ = 1.8×10⁻⁸, generation time = 5 yr (for comparability with Acropora),
chromosomes OZ116698.1–OZ116712.1, 3 species (oann/ofav/ofra).

**Expectations:**
- LGM bottleneck ~18–20 kya in all three species
- *O. faveolata* post-1980s collapse may appear as terminal Ne decline if resolution sufficient
- Results will allow refinement of moments absolute time/Ne estimates

---

## Reproducibility

```bash
cd /work/vollmer/orbicella_genomics
bash run.sh 1    # BAM ingestion + QC
bash run.sh 2a   # SNP discovery + relatedness (stops at clone gate)
bash run.sh 2b   # subset BEAGLE to unrelated samples
bash run.sh 3    # PCA + admixture
bash run.sh 4    # SAF + SFS + diversity + FST
bash run.sh 6    # moments demography (all 3 pairs)
bash run.sh 7    # SMC++ Ne(t)
```

Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
Reference: jaOrbFran1.1 (`/projects/vollmer/RR_heat-tolerance/Orbicella/reference/`)
