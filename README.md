# Strong three-almost-perfect matchings in trees

Exact verification code and recorded results for the research project **Trees maximizing strong three-almost-perfect matchings**, by Marek Gałązka and Hanna Wdowicka.

A matching pairs neighbouring vertices without using a vertex twice. Here the three uncovered vertices must be leaves of the original tree. The program counts these matchings, checks a weighted Gallai–Edmonds reduction, and compares an explicit five-group construction with generalized four-spiders.

This repository contains code, numerical outputs and reproduction instructions. The manuscript, its LaTeX sources and the cover letter are not included.

## Reproduce the checks

Tested with Python **3.12.14** and NetworkX **3.7**. Use Python 3.12 to reproduce the recorded environment. The pinned dependency requires Python 3.12 or later, excluding Python 3.14.1.

```bash
python -m venv .venv
# Activate the environment using the command appropriate for your shell.
python -m pip install -r requirements.txt
python verify.py --output run-results
```

Run Python without `-O` or `-OO`, because assertions implement the checks. The full run took about 10 seconds in the recorded environment; timing varies by machine. The output directory is created if needed and its result files are overwritten. `run-results` keeps the archived `results` directory intact.

For a smaller tree-enumeration run:

```bash
python verify.py --max-n 11 --output quick-results
```

This reduces the tree range, but still runs the allocation and explicit-family checks. Its summary should not be confused with the archived full run.

## What is checked

| Check | Archived count |
| --- | ---: |
| All nonisomorphic trees of orders 2 through 17 | 81,136 |
| Trees independently checked by deleting every leaf triple, through order 11 | 435 |
| Weighted-reduction and eligible-set instances, through order 12 | 9,628 |
| Nonnegative five-part allocations with total at most 40 | 1,221,759 |
| Explicit proportional constructions | 13 |
| Rounded optimizing constructions at six large orders | 7 |
| Failed assertions | 0 |

The independent leaf-triple count is a count of **trees tested**, not of individual triples. All permitted vertex sets are tested through order 7. At larger tested orders, the full, empty and original-leaf sets and five pseudorandom sets are used, with seed `20261030`. The reduction routines accept forests; the archived reduction tests use connected trees.

Three counting routes are implemented:

1. `coefficients`: rooted-forest polynomial dynamic programming. Coefficient `j` counts matchings leaving exactly `j` eligible vertices uncovered.
2. `direct_leaf_count`: delete each triple of original leaves and use NetworkX's general maximum-matching solver. A forest has at most one perfect matching.
3. `reduced_count`: construct Gallai–Edmonds sets from vertex-deletion matching numbers, compress the degree-two part, and enumerate branch-to-slot matchings directly. It does not use the specialized cubic formula.

All verification arithmetic is exact integer arithmetic.

## Output files

| File | Meaning |
| --- | --- |
| [`tree_enumeration.csv`](results/tree_enumeration.csv) | Exhaustive maxima, numbers of maximizers, and maxima restricted to deficiencies 1 and 3. |
| [`explicit_family.csv`](results/explicit_family.csv) | Counts for `D(2t,2t,t,2t,2t)` and the balanced positive-arm four-spider of the same order. |
| [`integer_formula.csv`](results/integer_formula.csv) | Values of the integer polynomial optimum `F(S)` and all maximizing middle weights. |
| [`summary.json`](results/summary.json) | Versions, check counts, elapsed time and SHA-256 of `verify.py`. |

The polynomial is

```text
P(a,b,y,c,d) = a*b*(c+d) + c*d*(a+b) + y*(a+b)*(c+d).
```

The construction `double_fork(a,b,y,c,d)` has a seven-vertex core and five groups of length-two pendant paths. For positive parameters its order is `2*(a+b+y+c+d)+7` and its count is `P`. Zero outer parameters can create extra original leaves, so this identity must not be applied to such parameters without rechecking the leaf condition.

In `integer_formula.csv`, `n=2*S+7`. The eventual extremal theorem applies only when `global_theorem_applies` is `True`, namely `n>=2007`. Rows below that threshold are polynomial optima, not asserted maxima over all trees. The threshold is sufficient and is not claimed to be optimal. A better count than a four-spider alone does not establish a global maximum.

In `tree_enumeration.csv`, a zero in a deficiency-specific column is also used when that class is absent. Even orders have count zero by parity. For order 3 there is no tree of deficiency 3.

Finite checks support the formulas and detect implementation mistakes. They do not prove the infinite extremal statements or establish novelty. The mathematical proofs belong to the separate manuscript.

## License and assistance

The code and recorded outputs are provided under the [MIT License](LICENSE).

The verification code was developed with assistance from ChatGPT (OpenAI).
