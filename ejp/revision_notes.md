# EJP revision for author review - 22 September 2026

EJP is the submission target.
No submission has been made by this workflow.

- Abstract and cover letter lead with logarithmic accuracy-horizon gains during exponential error growth. RK4 step and precision gains are qualified by the dominant error source.
- Conclusion proposes testing the same scaling with Lorenz or a driven pendulum; these are extensions, not reported experiments.
- Removed the unverified Calvao-Penna optimum-step comparison. Their citation remains for the general numerical-methods review described in the abstract.
- Shortened shadowing discussion and removed the unsupported guarantee of accurate long-time averages.
- Acknowledgments name Claude Code (Claude Opus 5 and Claude Fable 5.1, model IDs recorded in the prior session) and OpenAI Codex (GPT-6). The statement describes prose revision, computational-code assistance, and LaTeX work. Authors must review and accept responsibility for all content.
- Author order remains Arjun Pemmasani, Luke Wang.

## Validation

From ejp/: texcount -inc -sum manuscript.tex.
Text 3166; headings 74; captions/other float text 313. Total words: 3553,
including abstract and end matter, excluding bibliography. TeXcount's broader
sum is 3749, counting 191 inline and 5 displayed math expressions as units.
Shared numbers.tex is included; TeXcount is a source-word estimate, not a
count of every rendered mathematical symbol or generated table entry.
No fourth suggested problem exists in this EJP version, and no large cut was
needed to meet 4000 words.

Manuscript rebuilt with pdflatex, bibtex, pdflatex, pdflatex.
Cover letter rebuilt with pdflatex. Check final PDFs before upload.
No simulations were changed or rerun in this editorial revision.
Both authors should approve the final manuscript and disclosure before submission.

Final checks: manuscript 15 pages; cover letter 1 page. Rendered all manuscript pages and the letter and inspected their layout. No clipping observed. Final logs contain no unresolved-reference or overfull-box warnings; git diff --check passed. The revised files were prepared for author review.

## Follow-up corrections

Section 4.2 now gives gamma*8*ln(2)/lambda = 9.5060 using gamma=0.6
and lambda=0.35. Checked the data plotted by code/figures.py against
data/summary.json, precision_sweep, horizon["0.01"]. For p = 24, 28,
32, 36, 40, 48, 53, 64, 80, T(k=8)-T(k=16) is respectively
17.125, 8.8125, 3.875, 19.3125, 4.5, 12.8125, 6.6875, 18.875,
21.5625. Mean 12.6181; range 3.875-21.5625. The text reports approximate
agreement, not a rigid or exact offset. Figure 3(b) inspected.

Introduction names de Oliveira directly. Cover letter now explicitly states
both authors approved submission, as instructed by Arjun. Wording checked
against the author's Math131 homework 2 proofs and Math161 homework 1
solutions in Documents/LatexPSETS/latex-psets. The revision uses the same
relation-substitution-conclusion sequence. No results or figures changed.

Final follow-up validation: TeXcount 3216 text + 74 headings + 313 captions = 3603 words; broader sum 3803. Both PDFs rebuilt without unresolved-reference or overflow warnings. Manuscript 15 pages; cover letter 1 page. Supersedes the earlier counts above.
