# classification-qc

Quality checks for taxonomic classification output.

A metagenomic classification pipeline will happily return a confident-looking answer from a run that went wrong: host depletion that failed, a reference database missing the organism actually present, a reagent contaminant that swamped a low-biomass sample. The classifier does not know any of this. It reports what matched.

This is a small, dependency-free module that inspects classification output and flags the patterns that mean the result should not be trusted, before anyone acts on it.

## Install

No dependencies beyond the Python standard library (3.8+). Clone the repository and import it.

```bash
git clone https://github.com/<your-username>/classification-qc.git
cd classification-qc
python demo.py
```

## Use

```python
from classification_qc import grade_classification

reads = [
    {"read_id": "r1", "seq": "ACGT...", "tax_name": "Bacteroides fragilis"},
    {"read_id": "r2", "seq": "ACGT...", "tax_name": "Homo sapiens"},
    {"read_id": "r3", "seq": "ACGT...", "tax_name": "NA"},
]

for flag in grade_classification("gut", reads):
    print(flag.severity, flag.flag_type, flag.message)
```

Each read is a mapping with `read_id`, `seq` and `tax_name`, where `tax_name` is the assigned taxon or `"NA"` / `None` if the classifier assigned nothing.

## What it checks, and why

| Check | Severity | Scope | Rationale |
|---|---|---|---|
| `HOST_DNA_LEAK` | CRITICAL | per read | A read classified as *Homo sapiens* means host depletion did not fully work. This applies whatever the sample type, and it points at one specific upstream step. |
| `HIGH_UNCLASSIFIED_RATE` | CRITICAL | run | Above 50% unclassified, the problem is usually the pipeline rather than the biology: a reference database that does not cover what is present, or a poor run. |
| `UNEXPECTED_NICHE_TAXON` | WARNING | per read | A genus that does not fit the declared sample type is worth a look. Deliberately a warning, not an error: real biology is surprising, and the expectation set is small. |
| `DOMINANT_TAXON` | WARNING | run | Environmental samples are normally diverse. One organism above 80% of classified reads suggests a contamination spike or a processing artifact rather than community structure. |

Two of these are per-read and two are run-level, and the distinction is carried in the output: a `Flag` with `read_id=None` is a statement about the run, not about one observation. They call for different responses, so they should not look the same to whatever consumes them.

### Host DNA short-circuits the niche check

A human read is out of niche for every sample type, so the naive implementation reports both `HOST_DNA_LEAK` and `UNEXPECTED_NICHE_TAXON` for it. Both statements are true, but the second adds nothing and dilutes the first. The host check runs first and stops further evaluation of that read.

### An unknown sample type disables the niche check

`NICHE_EXPECTATIONS` covers `gut`, `skin` and `soil`. A sample type outside that mapping turns the niche check off rather than flagging every read in it. Flagging everything would be the same as flagging nothing, with more noise. Pass `niche_expectations=` to supply your own.

## Deliberate limitations

Worth stating plainly, because a QC tool that oversells itself is worse than none:

- **The niche expectation sets are illustrative, not authoritative.** Four genera per environment is a legible starting point for catching the obviously out-of-place, not a microbiome reference. Treat `UNEXPECTED_NICHE_TAXON` as a prompt to look, never as evidence of error.
- **Genus is taken as the first whitespace-separated token** of the taxon name. This holds for conventional binomials and breaks on strain designations and unconventional naming.
- **An empty flag list does not mean the classification is correct.** It means none of these four patterns were detected. These checks catch specific, common failure modes; they do not validate the assignment itself.
- **Thresholds are defaults, not findings.** 50% unclassified and 80% dominance are reasonable starting points, not values derived from a particular dataset. Both are parameters for that reason, and should be set against your own data.

## Possible extensions

Not implemented, and each would catch something the current checks miss:

- **Same sequence, different classification** within a run, which detects non-determinism or a pipeline bug rather than a sample problem.
- **Duplicate `read_id` with differing sequences**, an input-integrity check.
- **A dedicated contaminant blocklist** (PhiX, common reagent organisms) separate from the niche logic, since a known contaminant is a different claim from an implausible one.
- **Per-taxon read-count floors**, so single-read assignments are not treated with the same confidence as well-supported ones.

## Tests

```bash
python -m unittest discover -s tests -t . -v
```

19 tests covering each check, both severities, the per-read and run-level distinction, threshold boundaries, and the degenerate inputs: an empty read list, reads where every assignment is unclassified, and `NA` arriving as either the literal string or a null.

## Origin

The core of this was written in a timed one-hour exercise: four checks, a flag schema defined from scratch, and handling for the degenerate cases. This repository is that work made standalone, with configurable thresholds, the host-check short-circuit described above, a test suite, and the reasoning written down.

## License

MIT. See [LICENSE](LICENSE).
