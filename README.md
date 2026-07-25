# Split Iterated-Squaring Polynomials and Norm-like Gaussian Integers

Research on the **$E_n$ problem** — finding integers $(L_1,\dots,L_n)$ such that
the iterated-squaring polynomial $E_n(x)=((\cdots(x^2-L_1)^2-\cdots)^2-L_n$ splits
completely over $\mathbb{Z}$. An $E_5$ solution would give an ideal symmetric
Prouhet–Tarry–Escott solution of size 16, well beyond the current state of the art.

## Layout

- [`paper.md`](paper.md) — the paper (chapters 1–5 drafted; §§6–9 hold inline `TODO`s).
- [`TODO.md`](TODO.md) — workstreams (literature / theory / code / writing).
- [`drafts/sources_insight.md`](drafts/sources_insight.md) — synthesis of the source papers.
- [`drafts/verify_claims.py`](drafts/verify_claims.py) — reproduces every numeric claim in §§2–5.
- [`sources/`](sources/) — the reference PDFs.

## Build the PDF

The paper is Markdown with LaTeX math. Convert it with [pandoc](https://pandoc.org)
and a LaTeX engine:

```sh
pandoc paper.md -o paper.pdf \
  --pdf-engine=xelatex \
  -V geometry:margin=1in \
  -V linkcolor:blue
```

Notes on the flags:

- **`--pdf-engine=xelatex`** is required: the prose contains Unicode (`§`, `—`,
  `√`, `≈`, …) that the default `pdflatex` engine does not handle.
- **`-V geometry:margin=1in`** keeps every display formula inside the margins.
- **`-V linkcolor:blue`** colors the in-text citations, which are clickable links
  to the References section.

### Prerequisites

- `pandoc` (tested with 3.1.3)
- A TeX distribution providing `xelatex` (TeX Live: `sudo apt install texlive-xetex`,
  or the full `texlive-full`).

Sections are numbered manually in the source, so do **not** pass
`--number-sections` (it would double-number), and there is deliberately no
table of contents.

## Verify the computations

```sh
python3 drafts/verify_claims.py   # needs sympy
```

This checks the exemplar $E_4$ tower recovery and degree-16 polynomial identity,
the $\mathbb{Z}[\sqrt2]$ property of $L_{2,d}$, the divisibility ladder, the
$E_5$ lower bound $\sqrt{L_5}\ge C'_{16}/2$, and the §3 tower recursion.
