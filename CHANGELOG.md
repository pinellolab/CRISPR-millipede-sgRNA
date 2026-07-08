# Changelog

Entries here are not yet assigned a version — the user reviews accumulated changes and picks the next release number.

## [0.2.2] - 2026-07-08

### Performance

- **`get_mutation_profile`** per-guide linked-mutation encoding is now **vectorized**. The previous inner loop built a `pd.Series` per position and did a `pd.concat` per observed allele to prepend the `SequenceType` level, then a final `pd.concat(..., axis=1).transpose()` over all alleles; this is replaced by `get_substitution_encoding_batch`, which encodes all of a guide's aligned alleles at once into a single `M×6L` DataFrame via numpy broadcasting. **Bit-identical** to the previous output (validated: 729/729 protospacer-tier and 711/711 protospacer+surrogate+barcode-tier guides `.equals` on a real sample; `batch == per-allele` over 400 synthetic trials including `N`/`-`/`X` reference bases) and roughly **38× faster** (one sample's protospacer-tier profile: ~510 s → ~13 s).

### Added — backward-compatible keyword arguments (defaults reproduce previous behaviour)

- `compute_unlinked: bool = True` — when `False`, skip building the per-allele *unlinked* mutation tables (`determine_mutations_in_sequence`); the returned object's `all_observed_*_unlinked_mutations_df` stay `None`. Useful when only the `linked_mutations_whitelist_reporter_dict` is consumed; the linked dict is unchanged.
- `max_reads_per_guide: Optional[int] = None` / `max_alleles_per_guide: Optional[int] = None` — deterministic per-guide cap that keeps the most-abundant alleles first (cumulative reads ≥ N, or top-N alleles), to bound pathologically high-coverage guides. **Note:** the abundance cap biases per-variant (bystander-edit) frequencies — it is a bounded-work knob, not for frequency-sensitive use; leave unset (the default) for exact results.

## [0.2.1] - 2026-06-26

### Performance

- **`get_mutation_profile`** gains two backward-compatible keyword arguments (defaults reproduce the previous behaviour exactly):
  - `strategies: Optional[List[str]] = None` — restrict computation to the named ambiguity×UMI strategy base-names (e.g. `["ambiguous_accepted_umi_noncollapsed"]`); `None` computes all nine as before.
  - `cores: int = 1` — parallelise the per-guide loop across worker processes; `1` runs serially. The per-guide body is factored into a module-level worker (`_gmp_process_guide`) that reads the shared allele Series via fork copy-on-write; results are accumulated in input order so the output is identical to the serial path.
- On the equivalent change in the deployed older-API copy, computing a single strategy with `cores=12` was 713/713 guides bit-identical to the full nine-strategy serial output and cut one sample's mutation profile from >30 min to ~230 s.

## [0.2.0] - 2026-05-16

### Added — public API

- **`crispr_ambiguous_mapping.api`** with three stage-aligned entry points:
  - `map_fastq(library, fastq_r1_fns, fastq_r2_fns, *, config=ParsingConfig(...), **overrides)` — thin wrapper over the legacy `get_whitelist_reporter_counts_from_fastq`.
  - `count(result)` — returns the per-tier count-Series container.
  - `alleles(result, tier, *, contains_guide_surrogate, contains_guide_barcode, contains_guide_umi)` — wraps `get_matchset_alleleseries`; raises a clear `ValueError` on slim results.
- **`ParsingConfig`** dataclass — IDE-friendly bundle of the 50 parsing + threshold kwargs.
- **`MatchTier`** `str` enum (`PM`, `PM_SM`, `PM_BM`, `PM_SM_BM`, `PM_MISMATCH_SM`, `PM_MISMATCH_SM_BM`) — backward-compatible via `str` base class.
- **`save(result, directory)` / `load(directory)`** — parquet + JSON cross-language durable serialization (§7.4). Count series, QC summary, and `CountInput` round-trip; the per-observation inference dict stays pickle-only.
- **`crispr_correct`** top-level package alias — `import crispr_correct as cc` forward-looking name.

### Added — CLI

Flat-arg command-line entry point `crispr-correct` (§4.5). Every `ParsingConfig` field maps to `--flag value` (field names with `_` → `-`); no YAML, no config file.

- `crispr-correct map` — run mapping from FASTQs.
- `crispr-correct count` — emit one tier's count Series as TSV.
- `crispr-correct save` / `load` — round-trip a mapping result through a parquet/JSON directory.
- `crispr-correct alleles` — post-processing allele extraction to parquet.

### Added — validation + demo

- `tests/test_simulation_e2e.py` — 11-stage end-to-end sweep exercising every public surface (Python `map_fastq` / `count` / `alleles` / `save` / `load` + the five CLI subcommands) on the simulation fixtures with bit-equality and ground-truth assertions at every stage.
- `tests/simulation/phase8_validation_and_demo.ipynb` (wrapper folder) — 35-cell validation report + feature demo + guided tour of Phase 1-7 improvements with `git show --stat` snippets. Runs in ~90 s on the sim fixture. Outputs committed so it renders as a report on GitHub without running locally.
- **Fixed**: `MatchTier.__str__` was returning the Enum repr (`"MatchTier.PM_SM_BM"`) instead of the underlying string value, breaking `getattr(result, str(tier))`. Added `__str__` override returning `self.value`. Equality with plain strings via the `str` subclass is unaffected.
- **Fixed**: `crispr-correct alleles` built the strategy attribute name without the `ambiguous_` prefix, causing subprocess failures.

### Added — testing / CI

- GitHub Actions workflow `.github/workflows/ci.yml` runs smoke tests + 135-mode simulation regression (§7.7). Fast subset (~30 s) on push, full matrix on PR.
- `tests/fixtures/` checked in (library + 8k-read simulated FASTQs + truth parquet) so CI runs without external data.
- `tests/test_smoke.py` (13 tests) covers API surface, CLI, save/load round-trip.
- `tests/test_simulation.py` (135 parametrized comparisons across 8 parse modes × 6 tiers × 9 strategies).

### Changed

- **`MatchSetSingleInferenceMatchResultValue.matches`** is now `Optional[Tuple[Tuple[Any, ...], ...]]` (§2.4). Was `Optional[pd.DataFrame]` holding a per-observation iloc slice of the whitelist; now stores a plain tuple of positional row tuples. Per-match footprint drops ~0.8 KB → ~100 B; only matters when `retain_inference_results=True`. Column order matches whitelist registration (`protospacer` [, `surrogate`] [, `barcode`]). No compat shim — readers of `result.<tier>.matches` as a DataFrame must adapt.
- Padded whitelist DataFrame is released immediately after encoding (§2.7); was kept alive through the whole inference loop, duplicating the whitelist for the multiprocessing pool.
- **`retain_inference_results: bool = False` is the default** on `get_whitelist_reporter_counts_from_fastq` (§2.2). Result pickle shrinks ~15× (default) / ~93% at sim scale. Post-processing functions raise `ValueError` with an actionable message when called on a slim result.
- **`contains_surrogate` → `contains_guide_surrogate`** on `CountInput` and throughout post-processing kwargs (§4.3). Legacy `contains_surrogate` remains as a deprecated `@property` alias through the next release; removed after.
- **`contains_barcode` / `contains_umi` → `contains_guide_barcode` / `contains_guide_umi`** on `get_matchset_alleleseries`, `get_mutation_profile`, `tally_linked_mutation_count_per_sequence` (§1.4). No compat shim.
- Whitelist DataFrame columns `surrogate` / `barcode` are now accepted as `guide_surrogate` / `guide_barcode` (auto-renamed internally).
- `print()` statements replaced with module-level `logging.Logger` instances (§4.7 / §7.6). `logging.basicConfig(level=logging.INFO)` enables the verbose trace.
- Optional `tqdm` progress bar around the inference `pool.imap` (§4.8).
- Per-observation `pd.Series` construction replaced with a plain dict (§3.5).
- Redundant per-observation Hamming computations deduped — surrogate encoded once, full-whitelist Hamming computed once per read; subsets index-gather (§3.3).
- LUT-based DNA encoding replaces `np.vectorize` (§3.1).
- `pd.merge` in the per-observed-sequence loop replaced with `set.intersection` on tuple-of-guide-tuples (§3.4).
- O(N²) `.apply(axis=1)` cross-product in counter-series build replaced with `pd.DataFrame.from_records(counterdict.items())` — the counter-series stage went from ~21 317 s to ~53 s (~400×) on the AVITI 100k-read profile.
- `parse_fastq` collapsed from ~700 lines / 32 nested branches to ~150 lines / one generic component loop (§4.6).
- Default `surrogate_hamming_threshold_strict` corrected from 2 to 10.
- **Fixed** surrogate-length truncation bug (§1.1): observed surrogate is now clamped to the surrogate library length (32 bp) instead of the protospacer library length (20 bp). Output values change for any surrogate-involving tier — the scCRISPR golden pickle was re-baselined.

### Removed

- Deprecated `get_whitelist_reporter_counts_from_umitools_output` entry point and the two deprecated parsing modules (`reporter_umitools_fastq_parsing`, `guide_raw_fastq_parsing`) (§8). In-tree drivers must migrate to `get_whitelist_reporter_counts_from_fastq`.
- `non_error_dict` field on `QualityControlResult` tier objects (§2.1) — duplicated ~78% of the result pickle.
- `original_df` / `hamming_min_match_df` DataFrame fields on `HammingThresholdGuideCountError` subclasses (§2.3) — replaced with lightweight `n_whitelist_candidates` / `n_hamming_min_match` ints.
- `store_intermediates` flag — unused (§4.9).
- Legacy encoding helpers `encode_DNA_base_{whitelist,observed}{,_vectorized}`, `numpify_string{,_vectorized}` (§5.4) — superseded by the LUT-based `encode_DNA_sequence_*` / `encode_guide_series_*`.

### Performance deltas (cumulative Phase 1-5)

- **AVITI TCG** (100k reads, 1186-guide HBG library, full triplet + UMI): 272.3 s / 1461 MB → 155.3 s / 982 MB (**−43% wall, −33% peak RSS**).
- **chrX TGC** (100k reads, 11035-guide library, no UMI): 784.4 s / 2579 MB → 505.7 s / 1103 MB (**−36% wall, −57% peak RSS**).
- **scCRISPR pytest**: 75 s → 43 s (**−43% wall**).
- **Default result pickle**: 2.40 MB → 0.16 MB (**−93%**) at simulation scale.

### Dependencies

- Added `click ^8.1` (CLI).
- Added `pyarrow >=11,<22` (parquet save/load).
