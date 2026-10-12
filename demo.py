"""Three worked examples.

Run with:  python demo.py
"""

from classification_qc import grade_classification


def show(title, sample_type, reads):
    print(f"\n{title}")
    print(f"  sample_type={sample_type!r}, {len(reads)} reads")
    flags = grade_classification(sample_type, reads)
    if not flags:
        print("  no flags raised")
        return
    for f in flags:
        scope = "run" if f.is_run_level else f.read_id
        print(f"  [{f.severity:<8}] {f.flag_type:<24} ({scope}) {f.message}")


def read(read_id, tax):
    return {"read_id": read_id, "seq": "ACGT", "tax_name": tax}


# 1. A gut sample that looks the way a gut sample should.
show(
    "1. Plausible gut community",
    "gut",
    [read(f"a{i}", "Bacteroides fragilis") for i in range(4)]
    + [read(f"b{i}", "Bifidobacterium longum") for i in range(3)]
    + [read(f"c{i}", "Escherichia coli") for i in range(3)],
)

# 2. Host depletion did not fully work, and one read is out of niche.
show(
    "2. Host DNA surviving depletion",
    "gut",
    [read("h1", "Homo sapiens"), read("h2", "Homo sapiens")]
    + [read("x1", "Streptomyces coelicolor")]
    + [read(f"b{i}", "Bacteroides fragilis") for i in range(7)],
)

# 3. Low-biomass sample: most reads unclassified, and what did classify is
#    dominated by a single organism. Together these suggest the signal is
#    contamination rather than community structure.
show(
    "3. Low-biomass sample with a contamination spike",
    "skin",
    [read(f"u{i}", "NA") for i in range(7)]
    + [read(f"r{i}", "Ralstonia pickettii") for i in range(3)],
)

print()
