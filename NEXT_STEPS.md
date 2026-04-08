# Next Steps

Current state as of 2026-04-04.

---

## Pipeline status

| Segment | Status | Notes |
|---------|--------|-------|
| Seg 1: BAM ingestion + QC | ✅ Complete | 140 samples |
| Seg 2a: SNP discovery + relatedness | ✅ Complete | ~7.2M SNPs; 6 clones excluded; Ofra_26 mislabel detected |
| Seg 2b: Subset BEAGLE + grouped ngsRelate | ✅ Complete | 134 unrelated samples |
| Seg 3: LD + PCA + admixture | ✅ Complete | K=2 best (ΔK=26M); clean 3-species separation at K=3 |
| **Lineage assignment gate** | **⬜ Pending** | **Run lineage_assign.py before Seg 4** |
| Seg 4: Diversity + FST | ⬜ Pending | SAF, SFS, π, θ, Tajima's D; 3 inter-species FST |
| Seg 6: Demography (moments) | ⬜ Pending | All 3 inter-species pairs |
| Seg 7: SMC++ | ⬜ Pending | Ne(t) for all 3 species |

---

## Immediate next step

```bash
# On Discovery
cd /work/vollmer/orbicella_genomics
python /projects/vollmer/RRorbicella_angsd/workflow/scripts/lineage_assign.py
bash run.sh 4
```

---

## Key findings so far

- **K=2 best** by Evanno ΔK — primary signal is *O. annularis* vs *O. faveolata*+*O. franksi* clade
- **K=3** cleanly separates all 3 species; admixed individuals visible between *O. faveolata* and *O. franksi*
- **LD decays by ~8–10 kb** (r²<0.1) — notably faster than *Acropora*, suggesting larger Ne
- **Ofra_26** confirmed mislabel (*O. faveolata* genotype labeled as *O. franksi*) — notify data provider
