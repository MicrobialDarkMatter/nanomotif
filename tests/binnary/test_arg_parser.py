import pytest
from nanomotif.argparser import create_parser


def test_argparser_detect_contamination_with_specific_data(mocker):
    # Mock sys.argv with the specific command line you provided
    mocker.patch(
        "sys.argv",
        [
            "nanomotif",  # Assuming 'binnary.py' is a typo and should be 'binary.py' or similar
            "detect_contamination",
            "--motifs_scored",
            "data/motifs-scored.tsv",
            "--bin_motifs",
            "data/bin-motifs.tsv",
            "--contig_bins",
            "data/bins.tsv",
            "--out",
            "output",
        ],
    )

    # Parse the arguments
    parser = create_parser()
    args = parser.parse_args()

    # Assertions to verify that all arguments are correctly parsed
    assert args.command == "detect_contamination"
    assert args.motifs_scored == "data/motifs-scored.tsv"
    assert args.bin_motifs == "data/bin-motifs.tsv"
    assert args.contig_bins == "data/bins.tsv"
    assert args.out == "output"
