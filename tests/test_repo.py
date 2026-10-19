"""Static repo checks (stdlib only; no datasets or scientific env required).

Run:  python -m unittest discover tests
"""
import ast
import glob
import json
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NOTEBOOKS = sorted(glob.glob(os.path.join(ROOT, "*", "*.ipynb")))


def code_cells(path):
    with open(path) as fh:
        nb = json.load(fh)
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


class NotebookTests(unittest.TestCase):
    def test_found_notebooks(self):
        self.assertGreater(len(NOTEBOOKS), 20)

    def test_code_cells_parse(self):
        for nb in NOTEBOOKS:
            for i, src in enumerate(code_cells(nb)):
                lines = [l for l in src.split("\n") if not l.lstrip().startswith(("%", "!"))]
                try:
                    ast.parse("\n".join(lines))
                except SyntaxError as e:
                    self.fail(f"{nb} cell {i}: {e}")

    def test_data_dir_is_repo_root(self):
        """Notebooks live one level down; data must go to ../data (see clean.sh)."""
        for nb in NOTEBOOKS:
            for src in code_cells(nb):
                for m in re.finditer(r'^DATA_DIR\s*=\s*["\']([^"\']+)["\']', src, re.M):
                    self.assertEqual(m.group(1), "../data", nb)
                for m in re.finditer(r'["\'](\.?/?data/[^"\']*)["\']', src):
                    self.fail(f"{nb}: bare relative data path {m.group(1)!r}")

    def test_savefig_dirs_created(self):
        for nb in NOTEBOOKS:
            src = "\n".join(code_cells(nb))
            if "FIG_DIR" in src and "savefig" in src:
                self.assertRegex(src, r"makedirs\(FIG_DIR", nb)


class DocTests(unittest.TestCase):
    def test_markdown_relative_links_resolve(self):
        for md in glob.glob(os.path.join(ROOT, "**", "*.md"), recursive=True):
            if "node_modules" in md:
                continue
            with open(md) as fh:
                text = fh.read()
            for m in re.finditer(r"\]\(([^)#\s]+)(?:#[^)]*)?\)", text):
                target = m.group(1)
                if re.match(r"[a-z]+:", target):
                    continue
                self.assertTrue(
                    os.path.exists(os.path.join(os.path.dirname(md), target)),
                    f"{md}: broken link {target}",
                )

    def test_clean_sh_covers_perturb_seq_outputs(self):
        clean = open(os.path.join(ROOT, "perturb_seq", "clean.sh")).read()
        for nb in glob.glob(os.path.join(ROOT, "perturb_seq", "*.ipynb")):
            src = "\n".join(code_cells(nb))
            for name in re.findall(r'"(\w+\.(?:h5ad|pkl|csv))"', src):
                self.assertIn(name, clean, f"{nb} writes/reads {name} not in clean.sh")


if __name__ == "__main__":
    unittest.main()
