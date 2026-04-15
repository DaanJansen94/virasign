# Z-score background correction

Virasign can compute a **background-corrected Z-score** per reported virus using **water controls** (e.g. negative controls). This mirrors a common idea in metagenomic reporting: quantify whether a signal is unusually high compared to background contamination.

This concept is used in CZ ID / IDseq as part of their “background model” reporting. See:
- IDseq paper (GigaScience, 2020): `https://ncbi.nlm.nih.gov/pmc/articles/PMC7566497/` (background models and z-scores described at a high level)
- CZ ID workflows wiki: `https://github.com/chanzuckerberg/czid-workflows/wiki`

---

## Why use a Z-score?

Some taxa show up at low levels in many runs due to:
- reagent / extraction contamination
- index hopping / low-level carryover
- laboratory environment background

By comparing each virus to your **water controls in the same run**, the Z-score answers:

> “Is this virus higher than what we typically see in water controls?”

---

## When Virasign computes it

- Z-scores are computed only when **≥2 water controls** are available.
- Water controls are auto-detected when the sample name contains `water`, `h2o`, or `h20` (case-insensitive).
- You can override auto-detection with **exact input FASTQ paths**:

```bash
virasign -i input_dir --zscore-controls /path/to/water1.fastq.gz,/path/to/water2.fastq.gz
```

When `--zscore-controls` is provided, Virasign **does not** use auto-detection (so you can exclude some “H2O_*” samples intentionally).

---

## What signal is used

Virasign computes the Z-score using the per-hit **remapped** `mapped_reads` (the same value reported in the final per-sample JSON).

To make the statistic stable for sparse count data, Virasign uses:

```text
y = log10(mapped_reads + 1)
```

---

## Formula

For a given virus (practical label) \(v\), and a set of water controls \(W\):

- Compute the transformed values in each water control:
  \[
  y_{w,v} = \log_{10}(x_{w,v} + 1)
  \]
  where \(x_{w,v}\) is `mapped_reads` for virus \(v\) in water sample \(w\). If a virus is absent from a water sample, \(x_{w,v}=0\).

- Compute background mean and standard deviation:
  \[
  \mu_v = \text{mean}(\{y_{w,v}\})
  \]
  \[
  \sigma_v = \text{stdev}(\{y_{w,v}\}) \quad (\text{sample stdev, ddof}=1)
  \]

For a non-water sample \(s\):

\[
z_{s,v} = \frac{y_{s,v} - \mu_v}{\sigma_v}
\]

### Zero-variance case (\(\sigma_v = 0\))

If all controls have exactly the same value for a virus, the classical Z-score is undefined. Virasign keeps the output numeric and directional by using a tiny \(\epsilon\) in the denominator and capping extremes:

- equal to controls \(\Rightarrow z=0\)
- higher than controls \(\Rightarrow\) large positive Z
- lower than controls \(\Rightarrow\) large negative Z

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
- `results_summary_*.html` and `results_summary_*.csv` include a **Z-score column**.

