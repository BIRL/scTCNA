"""scTCNA: single-cell transcriptional noise atlas for cancer.

A single-cell RNA-seq counterpart of The Cancer Noise Atlas (TCNA; Shah et al.,
J. Mol. Biol. 2026). Metric definitions follow the TCNA manuscript where it
defines them (CV, CV^2, mean absolute deviation, SD, differential noise as the
log2 ratio of CVs, DEPTH and DEPTH2); single-cell specific choices are
documented in METHODS.md and marked "scTCNA assumption" in the code.
"""

__version__ = "0.1.1"
