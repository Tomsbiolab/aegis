"""Tests for aegis orthology CLI helpers and edge cases."""

import pandas as pd

from aegis.cli.orthology import get_tiered_cardinality


class TestGetTieredCardinality:
    def test_basic_cardinalities_and_name_order_symmetry(self):
        # 1:1 (x1-y1), 1:N (x2-y2, x2-y3), N:1 (x3-y4, x4-y4), N:N path (x5-y5, x5-y6, x6-y6)
        rows = [
            ("x1", "y1"),
            ("x2", "y2"), ("x2", "y3"),
            ("x3", "y4"), ("x4", "y4"),
            ("x5", "y5"), ("x5", "y6"), ("x6", "y6"),
        ]
        expected = ["1:1", "1:N", "1:N", "N:1", "N:1", "N:N", "N:N", "N:N"]

        for name_a, name_b in (("Alpha", "Beta"), ("Zeta", "Beta")):
            df = pd.DataFrame({
                "annotation_A": [name_a] * len(rows),
                "annotation_B": [name_b] * len(rows),
                "gene_id_A": [a for a, _ in rows],
                "gene_id_B": [b for _, b in rows],
                "summary_score": ["high_confidence"] * len(rows),
            })
            assert list(get_tiered_cardinality(df, ["high_confidence"])) == expected

    def test_reversed_direction_flips_1N_and_N1_symmetrically(self):
        rows = [
            ("x1", "y1"),
            ("x2", "y2"), ("x2", "y3"),
            ("x3", "y4"), ("x4", "y4"),
            ("x5", "y5"), ("x5", "y6"), ("x6", "y6"),
        ]
        df_fwd = pd.DataFrame({
            "annotation_A": ["Alpha"] * len(rows),
            "annotation_B": ["Beta"] * len(rows),
            "gene_id_A": [a for a, _ in rows],
            "gene_id_B": [b for _, b in rows],
            "summary_score": ["high_confidence"] * len(rows),
        })
        df_rev = pd.DataFrame({
            "annotation_A": ["Beta"] * len(rows),
            "annotation_B": ["Alpha"] * len(rows),
            "gene_id_A": [b for _, b in rows],
            "gene_id_B": [a for a, _ in rows],
            "summary_score": ["high_confidence"] * len(rows),
        })
        both = pd.concat([df_fwd, df_rev], ignore_index=True)
        labels = list(get_tiered_cardinality(both, ["high_confidence"]))
        fwd_labels = labels[:len(rows)]
        rev_labels = labels[len(rows):]
        flip = {"1:1": "1:1", "1:N": "N:1", "N:1": "1:N", "N:N": "N:N"}
        assert [flip[x] for x in fwd_labels] == rev_labels
