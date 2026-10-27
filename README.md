Écriture

# Propagation of the Evoked Response Depending on Gamma State (Low / High)

This repository contains the notebooks extending the paper _“Responsiveness of cerebral cortex to hippocampal inputs depends on brain region and brain state.”_ The top of each notebook reproduces the original paper's notebook unchanged. All additions are **at the bottom**.

| Notebook | Contents |
| --- | --- |
| `00_ARTICLE.ipynb` | The original paper's notebook, unchanged. |
| `01_ARTICLE_BIDIR.ipynb` | A 10-node mean-field chain, **feedforward and bidirectional**, with 8 values of C, stimulation at node 0 or 5, followed by a **matched control** (same simulation, same noise, no stimulus). |
| `02_TVB_ARTICLE.ipynb` | The same setup implemented using **TVB** (the lab's software stack), plus the lab's **2D grid**. |
| `03_HUGO_BIDIR.ipynb` | Same as 01, using **Hugo's transfer function**. |
| `04_TVB_HUGO.ipynb` | Same as 02, using Hugo's transfer function. |

The notebooks are delivered executed: figures are visible without rerunning anything. They are also in figures_gamma/00 to figures_gamma/04, one folder per notebook.

## What Is Measured

- **Stimulus and measurement: as in the paper.** 200 Gaussian pulses (σ = 10 ms, 0.9 Hz).
- **Matched control.** Each condition is simulated twice with exactly the same noise: once with the stimulus and once without it.

## Figures at the Bottom of Each Notebook

1. Peak response at each node for each C value, in the format used in the paper.
2. Stimulus − Control peak response, expressed as a percentage of the control baseline (bars show the standard error across 200 trials).
3. Dependence on C at distances of 0, 2, and 4 from the stimulated node.
4. Time × distance plot.
5. Low and High traces with their corresponding controls, at distances of 0, 2, and 4, with the peak value displayed.
6. A final table showing the stimulus and control peak responses at distances of 0, 2, and 4.
