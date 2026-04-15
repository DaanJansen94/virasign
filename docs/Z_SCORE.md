# Z-score background correction

Virasign can compute a background-corrected Z-score per reported virus using water controls (e.g. negative controls). This mirrors a common idea in metagenomic reporting: quantify whether a signal is unusually high compared to background contamination.

---

## Why use a Z-score?

Some taxa show up at low levels in many runs due to:
- reagent / extraction contamination
- index hopping / low-level carryover
- laboratory environment background

By comparing each virus to your water controls in the same run, the Z-score answers:

> “Is this virus higher than what we typically see in water controls?”

---

## When Virasign computes it

- Z-scores are computed only when ≥2 water controls are available.
- Water controls are auto-detected when the sample name contains `water`, `h2o`, or `h20` (case-insensitive).
- You can override auto-detection with exact input FASTQ paths:

```bash
virasign -i input_dir --zscore-controls /path/to/water1.fastq.gz,/path/to/water2.fastq.gz
```

When `--zscore-controls` is provided, Virasign does not use auto-detection (so you can exclude some “H2O_*” samples intentionally).

---

## What signal is used

Virasign computes the Z-score using the per-hit remapped `mapped_reads` (the same value reported in the final per-sample JSON).

---

## Formula

Virasign computes, for each virus, a Z-score as the number of standard deviations that the sample’s log-transformed mapped_reads signal is above/below the mean of the log-transformed mapped_reads values in the selected water controls.

This follows the same background-correction idea used by CZ ID / IDseq background models. For more information, see:

- IDseq paper (GigaScience, 2020): `https://ncbi.nlm.nih.gov/pmc/articles/PMC7566497/`
- CZ ID workflows wiki: `https://github.com/chanzuckerberg/czid-workflows/wiki`

---

## How to interpret Z-scores

Z-scores are “standard deviation units” above/below background:

- **Z ≈ 0**: similar to water background
- **Z < 0**: lower than the water background mean
- **Z > 0**: above water background
- **Z ≥ 2**: ≥2 standard deviations above the control mean.
- **Z ≥ 3**: ≥3 standard deviations above the control mean.
- Very large Z (e.g. **50–100**) can happen when controls are very consistent (low variance) and the sample is higher than background.

---

## Where it appears in outputs

- Per-hit fields in `*_final_selected_references.json`:
  - `zscore`: the computed value
  - `zscore_controls`: list of water samples used for the background model
- `results_summary_*.html` and `results_summary_*.csv` include a Z-score column.

