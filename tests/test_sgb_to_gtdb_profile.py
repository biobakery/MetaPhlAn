import pytest

from metaphlan.utils.sgb_to_gtdb_profile import get_gtdb_profile, normalize_abundances


def test_normalize_abundances_merged_without_unclassified():
    abundances = {
        "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus reuteri": [
            60.0,
            25.0,
        ],
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis": [
            20.0,
            75.0,
        ],
    }

    normalized = normalize_abundances(abundances, merged_profiles=True, unclassified_fraction=0)

    values = list(normalized.values())
    assert values[0] == pytest.approx([75.0, 25.0])
    assert values[1] == pytest.approx([25.0, 75.0])


def test_get_gtdb_profile_merged_without_unclassified(tmp_path):
    mapping = tmp_path / "mpa_test_SGB2GTDB_r226.tsv"
    mapping.write_text(
        "SGB1\td__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus reuteri\n"
        "SGB2\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis\n"
    )

    merged_input = tmp_path / "merged_profile.tsv"
    merged_input.write_text(
        "#mpa_test\n"
        "clade_name\tsample1\tsample2\n"
        "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Limosilactobacillus|s__Limosilactobacillus_reuteri|t__SGB1\t60.0\t25.0\n"
        "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB2\t20.0\t75.0\n"
    )

    gtdb_output = tmp_path / "gtdb_merged.tsv"
    get_gtdb_profile(
        str(merged_input),
        str(gtdb_output),
        merged_profiles=True,
        gtdb_assignment_file=str(mapping),
    )

    rows = {}
    for line in gtdb_output.read_text().splitlines():
        if line.startswith("#") or line.startswith("clade_name"):
            continue
        parts = line.split("\t")
        rows[parts[0]] = [float(x) for x in parts[1:]]

    assert rows["d__Bacteria"] == pytest.approx([100.0, 100.0])
    assert rows[
        "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus reuteri"
    ] == pytest.approx([75.0, 25.0])
    assert rows[
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis"
    ] == pytest.approx([25.0, 75.0])


def test_get_gtdb_profile_merged_with_unclassified(tmp_path):
    mapping = tmp_path / "mpa_test_SGB2GTDB_r226.tsv"
    mapping.write_text(
        "SGB1\td__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus reuteri\n"
    )

    merged_input = tmp_path / "merged_profile_unclass.tsv"
    merged_input.write_text(
        "#mpa_test\n"
        "clade_name\tsample1\tsample2\n"
        "UNCLASSIFIED\t20.0\t50.0\n"
        "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Limosilactobacillus|s__Limosilactobacillus_reuteri|t__SGB1\t80.0\t50.0\n"
    )

    gtdb_output = tmp_path / "gtdb_merged_unclass.tsv"
    get_gtdb_profile(
        str(merged_input),
        str(gtdb_output),
        merged_profiles=True,
        gtdb_assignment_file=str(mapping),
    )

    rows = {}
    for line in gtdb_output.read_text().splitlines():
        if line.startswith("#") or line.startswith("clade_name"):
            continue
        parts = line.split("\t")
        rows[parts[0]] = [float(x) for x in parts[1:]]

    assert rows["UNCLASSIFIED"] == pytest.approx([20.0, 50.0])
    assert rows["d__Bacteria"] == pytest.approx([80.0, 50.0])


def test_get_gtdb_profile_single_without_unclassified(tmp_path):
    mapping = tmp_path / "mpa_test_SGB2GTDB_r226.tsv"
    mapping.write_text(
        "SGB1\td__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus reuteri\n"
    )

    single_input = tmp_path / "single_profile.tsv"
    single_input.write_text(
        "#mpa_test\n"
        "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
        "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Limosilactobacillus|s__Limosilactobacillus_reuteri|t__SGB1\t2|1239|91061|186826|33958|2742598|1598|\t80.0\t\n"
    )

    gtdb_output = tmp_path / "gtdb_single.tsv"
    get_gtdb_profile(
        str(single_input),
        str(gtdb_output),
        merged_profiles=False,
        gtdb_assignment_file=str(mapping),
    )

    rows = {}
    for line in gtdb_output.read_text().splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        rows[parts[0]] = float(parts[1])

    assert rows["d__Bacteria"] == pytest.approx(100.0)
