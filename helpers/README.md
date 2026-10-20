# helpers/

Claude-authored utility code for the notebooks, kept separate from `src/`
(which Nobu writes and maintains himself). Nothing here replaces or
duplicates `src/` logic -- it only wraps/combines existing `src/` functions
for a specific notebook step. If something here turns out to be broadly
useful, it should graduate into `src/` deliberately, not by accretion.

Currently empty: Step 1 of the October 2026 Kirishima rebuild
(notebooks/Elastic2DtopoBenchmark_oct2026.ipynb) only calls existing
src/ functions directly, so nothing has been needed here yet.
