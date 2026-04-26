# Pluggable Annotators + VEP REST (v0.12)

- **Status**: Accepted, executing
- **Targets**: v0.12.0-pilot
- **Roadmap item**: walkthrough §8 #1 (Real VEP for phase 1) reinterpreted under user-stated principles (API-first, modular, disease-agnostic, no 50 GB local caches)

## 1. Goal

Replace the curated `GENE` + `CONSEQUENCE` shortcut in phase 1 with a real Ensembl VEP call, via the public REST API rather than a local 50 GB cache. Do it through a pluggable annotator interface so future backends (OncoKB, ClinVar, AlphaMissense, local VEP CLI) slot in the same way.

## 2. Non-goals

- Local VEP CLI install. Defer until a real motivation surfaces (rate-limit pain, plugin requirements, regulatory reproducibility).
- OncoKB layer. Separate v0.13 work.
- Tier classification in phase 5. Separate v0.14 work.
- Mutation calling (Mutect2, Strelka). Separate; we still take pre-called VCFs as input.
- Pathogenicity score axis (REVEL, CADD, AlphaMissense). Future annotator backends; not in this release.

## 3. Constraints

- Stdlib only (`urllib.request`, `json`). No new dependencies.
- Em-dash gate must stay clean.
- Scorecard 10/10 PASS must hold throughout (toys default to curated path; their CONSEQUENCE values stay hand-baked).
- CI must not depend on network: VEP REST tests use mocked HTTP.
- Real-cohort runs (`bin/run_target_rt.py`) use VEP REST by default.

## 4. Architecture

```
bin/phase1_annotate.py
        |
        | --annotator curated|vep_rest  (default: curated)
        v
bin/annotators/
    base.py          Annotator protocol + Variant/VariantAnnotation dataclasses + factory
    curated_vcf.py   Reads gene/consequence/hgvsp_short from VCF INFO (current behavior)
    vep_rest.py      Batched POST to rest.ensembl.org/vep/human/region; on-disk cache
```

Phase 1 unchanged in spirit: it reads the VCF, calls the annotator to fill in `(gene, consequence, hgvsp_short)`, then runs the existing `classify_snv` rule. The annotator pattern is the only new abstraction.

## 5. Annotator protocol

```python
@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    info: dict[str, str]   # passthrough of VCF INFO key/value (curated backend reads from this)


@dataclass(frozen=True)
class VariantAnnotation:
    gene: str
    consequence: str       # Sequence Ontology term (e.g., "missense_variant")
    hgvsp_short: str       # Short form (e.g., "V550L"); empty string if not applicable
    source: str            # Annotator name that produced this row (e.g., "vep_rest")


class Annotator(Protocol):
    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        """One annotation per input variant, same order."""


def get_annotator(name: str, **opts) -> Annotator:
    """Factory. Accepted names: 'curated', 'vep_rest'."""
```

Per-variant degradation: if a backend cannot produce a real annotation for a given variant, it returns a `VariantAnnotation` with empty `consequence` (and best-effort `gene`/`hgvsp_short`). Phase 1's existing `classify_snv` already handles empty consequence by routing to PASSENGER, so no new failure paths.

## 6. CuratedVCFAnnotator

Wraps current behavior:

- `gene` from `INFO["GENE"]`
- `consequence` from `INFO["CONSEQUENCE"]`
- `hgvsp_short` from regex on `INFO["NOTE"]` (the existing `extract_hotspot`)
- `source = "curated"`

This is a literal extraction of the existing inline code in `bin/phase1_annotate.py`. The whole point is bit-for-bit equivalence on toy fixtures.

## 7. VEPRestAnnotator

- Endpoint: `POST https://rest.ensembl.org/vep/human/region`
- Headers: `Content-Type: application/json`, `Accept: application/json`
- Request body: `{"variants": ["<chrom> <pos> . <ref> <alt> . . ."]}` per the Ensembl VCF-style format. Batches of up to 200 variants per request (Ensembl's documented limit).
- Response: list of objects, each with a `transcript_consequences` array. Pick the canonical entry (first one with `canonical: 1`, fall back to the most severe `consequence_terms` when no canonical flag is set).
- HGVSp comes from `transcript_consequences[i].hgvsp`, format `ENST...:p.Val550Leu`. Helper converts to `V550L` (strip transcript prefix, drop `p.`, run 3-letter-to-1-letter table).
- `consequence`: first element of `consequence_terms`.
- `gene`: `transcript_consequences[i].gene_symbol`.
- `source = "vep_rest"`.

Cache:

- Default location: `data/vep_cache/<chrom>_<pos>_<ref>_<alt>.json`. Override via `--vep-cache-dir`.
- Each cached file is the raw VEP response object for that single variant (parsed from the batched response).
- Cache hit -> no API call. Cache miss -> request batched along with other misses, then each result written to its own cache file.
- `data/` is gitignored (already), so cache stays local; refresh is just `rm -rf data/vep_cache`.

Failure modes:

- HTTP 4xx/5xx from VEP -> log a warning to stderr, return a degraded `VariantAnnotation` (empty consequence).
- Network unreachable -> log warning, degraded annotations for the whole batch.
- Per-variant parse failure (e.g., variant on a non-canonical chromosome) -> degraded annotation for that one variant only.

Rate limit: Ensembl publicly enforces 15 req/sec, 55,000 req/hour. With 200 variants per batch, that's effectively unbounded for any cohort we'd care about (36 samples ~10 vars each = 1 batch). No in-process rate limiting needed yet.

## 8. Phase 1 wiring

`bin/phase1_annotate.py` adds:

- `--annotator {curated, vep_rest}` (default: `curated`).
- `--vep-cache-dir <path>` (default: `<repo>/data/vep_cache`).
- `--disease <name>` (default: `RMS`). Currently unused; reserved for future OncoKB/AlphaMissense backends that condition on disease. Stored in the annotator factory call so backends can stash it for their own use.

The existing `annotate_vcf()` function refactors to:

1. Parse VCF lines into `Variant` objects (chrom/pos/ref/alt/info).
2. Call `annotator.annotate_batch(variants)`.
3. For each `(variant, annotation)` pair, run `classify_snv(annotation.consequence, gene_kb, annotation.hgvsp_short)` exactly as today, populate the row dict.

Output column contract is unchanged; downstream phases are not affected.

## 9. Real-cohort wiring

`bin/run_target_rt.py` adds `"--annotator", "vep_rest"` to its phase 1 invocation. Net effect on a fresh cohort run: ~1 batched VEP call upfront (one per study fetcher run; maybe a few if variants are spread across many regions). Subsequent runs hit the cache and are instant.

## 10. Tests

`tests/test_annotators.py`:

- `TestCuratedVCFAnnotator`: regression. Uses an in-memory `Variant` with hand-built INFO dict; asserts the returned `VariantAnnotation` matches the inline parsing in `bin/phase1_annotate.py:annotate_vcf`.
- `TestVEPRestAnnotator`: monkey-patches `urllib.request.urlopen` to return a canned VEP response (loaded from `tests/data/vep_rest_responses/<variant>.json`). Asserts:
    - HGVSp 3-letter-to-1-letter conversion (e.g., `ENST00000379368.4:p.Val550Leu` -> `V550L`)
    - Canonical transcript selection
    - Multi-variant batch order preservation
    - Cache hit on second call (no HTTP invocation)
    - Cache write on first call
- `TestFactory`: `get_annotator("curated")` and `get_annotator("vep_rest")` return the right types; unknown names raise.

Phase 1 regression: existing toy fixtures run through `classify_snv` -> identical output. Verified by running the scorecard before and after the refactor.

## 11. CI

Adds a unittest invocation alongside the existing `test_cohort_visualize` step:

```yaml
- name: Run annotator unit tests
  run: |
    python3 -m unittest tests.test_annotators -v
```

CI does not hit the network because `TestVEPRestAnnotator` mocks urlopen.

## 12. Walkthrough + ADR

- `docs/adr/0003-pluggable-annotators.md`: records the abstraction and the decision to default to API-first backends with on-disk caching rather than 50 GB reference downloads.
- `docs/walkthrough.md` §2.1 and §8 #1 update to reflect the v0.12 reality.

## 13. Risks

- **VEP REST returns slightly different consequence terms than what `classify_snv` accepts.** Already mitigated: `PROTEIN_ALTERING_CONSEQUENCES` includes both forms (`missense` and `missense_variant`).
- **VEP HGVSp format edge cases**: insertions, deletions, frameshifts have non-`Val550Leu` formats. Mitigation: helper returns empty `hgvsp_short` for anything that doesn't match the simple substitution pattern. Phase 1's hotspot match silently fails to match (which is the right behavior).
- **Cache key collision with multi-allelic sites.** Mitigation: cache key includes all of (chrom, pos, ref, alt); two ALT alleles at the same site write to different files.
- **API-deprecation risk.** Mitigated by the abstraction itself: if Ensembl REST goes away, drop in a different backend without touching phase 1.

## 14. Acceptance

- `python3 -m unittest tests.test_annotators -v` exits 0.
- `python3 -m unittest tests.test_cohort_visualize -v` exits 0 (regression).
- `python3 bin/check_case_studies.py --quiet` exits 0 (scientific gate untouched).
- `python3 bin/run_target_rt.py` runs through, phase 1 calls VEP REST, cache populates, real cohort markdown contains hgvsp_short values for FGFR4 V550L etc.
- Tracked-markdown em-dash gate clean.
