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

---

## Relatedness and Clonality

133 unrelated samples retained after clone/relatedness filtering. Key findings:
- Ofra_26 is a clone of Ofav_106 (KING=0.493); excluded
- Ofra_6 is a putative *franksi* × *faveolata* hybrid (8.4% *faveolata* at K=3); excluded

---

## Population Structure

K=3 cleanly separates all three species. ΔK supports K=3 as the optimal clustering.
PCAngsd and NGSAdmix results are concordant. Lineage assignments at Q > 0.80 threshold used
for all downstream analyses (demography, FST).

**Note:** The genetic distance / UPGMA tree places *O. annularis* closer to *O. franksi*,
which is a reference genome bias artefact (alignment to *O. franksi* reference inflates
apparent similarity). Fukami et al. (2004) phylogeny (*O. faveolata* + *O. franksi* as
sisters) is supported by the demographic inference (see below).

---

## Genetic Diversity

| Species | π (genome-wide) | Tajima's D |
|---------|----------------|------------|
| *O. annularis* | 0.01241 | −0.59 |
| *O. faveolata* | 0.01228 | −0.20 |
| *O. franksi* | 0.01272 | −0.52 |

All three species have similar diversity (π ≈ 0.012–0.013). All Tajima's D values are
negative, consistent with recent population expansion or purifying selection. *O. faveolata*
has the least negative D (closest to 0), suggesting a more stable recent history relative
to the other two species.

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
*O. franksi* converged with T2 ≈ 0 (i.e., ancient migration collapsed to zero),
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
Oann→Ofav migration is ~3× stronger than Ofav→Oann (scaled m: 1.38 vs 0.47). This
pattern may reflect that *O. annularis* and *O. faveolata* have overlapping ranges and
spawning windows sufficient to maintain gene flow, even as they remain genetically
distinct species (FST pending).

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
```

Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
Reference: jaOrbFran1.1 (`/projects/vollmer/RR_heat-tolerance/Orbicella/reference/`)
