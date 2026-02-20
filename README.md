# Sequence Alignment

A Python implementation of optimal DNA sequence alignment using dynamic programming, with both a standard full-table approach and a memory-efficient linear-space variant.

---

## Problem

Given two DNA sequences **X** and **Y**, find the optimal alignment that minimises the total cost of:
- **Gap insertions** (`–`) — each gap incurs a fixed penalty **δ**
- **Mismatches** — each mismatched character pair incurs a substitution cost **α(x, y)**

An alignment maps each character in X and Y to either a matching position or a gap, for example:

```
X: A C G – T A
Y: A – G C T A
       ↑   ↑
     gaps  match
```

This is the classic **global sequence alignment** problem (Needleman–Wunsch), solved here in two ways.

---

## Algorithms

### 1. Basic Dynamic Programming (`basic_3.py` / `basic.sh`)

The standard **O(mn)** time and **O(mn)** space Needleman–Wunsch algorithm.

- Fills a full `(m+1) × (n+1)` cost table.
- Recovers the optimal alignment by backtracking through the table.
- Simple and easy to follow, but memory usage grows quadratically with sequence length.

**Recurrence:**

```
OPT(i, 0) = i · δ
OPT(0, j) = j · δ
OPT(i, j) = min(
    α(xᵢ, yⱼ)  + OPT(i-1, j-1),   # match / mismatch
    δ           + OPT(i-1, j),      # gap in Y
    δ           + OPT(i,   j-1)     # gap in X
)
```

### 2. Memory-Efficient DP (`efficient_3.py` / `efficient.sh`)

A **Hirschberg-style** divide-and-conquer algorithm that achieves **O(mn)** time with only **O(m + n)** space.

- Uses the observation that the optimal alignment score can be computed using only two rows at a time.
- Recursively splits the problem at the midpoint of Y and finds the optimal split point in X using forward and backward passes.
- Reconstructs the full alignment without ever storing the complete DP table.

| | Basic DP | Efficient DP |
|---|---|---|
| **Time** | O(mn) | O(mn) |
| **Space** | O(mn) | O(m + n) |
| **Best for** | Short sequences, clarity | Long sequences |

---

## File Structure

```
.
├── basic_3.py       # Basic O(mn) space DP implementation
├── basic.sh         # Shell script to run basic DP with example inputs
├── efficient_3.py   # Memory-efficient Hirschberg-style implementation
└── efficient.sh     # Shell script to run efficient DP with example inputs
```

---

## Usage

### Run directly with Python

```bash
# Basic DP
python basic_3.py <sequence_X> <sequence_Y>

# Memory-efficient DP
python efficient_3.py <sequence_X> <sequence_Y>
```

**Example:**
```bash
python basic_3.py ACGTAC ACGTTC
```

**Output:**
```
Alignment cost: 30
A C G T A C
A C G T T C
Time: 0.0012s  Memory: 0.18 MB
```

### Run via shell scripts

```bash
bash basic.sh
bash efficient.sh
```

The shell scripts run a suite of test cases and print the alignment cost, aligned sequences, execution time, and memory usage.

---

## Input Format

Sequences consist of the four DNA bases: `A`, `C`, `G`, `T`.

Costs are defined inside the Python files:

```python
# Gap penalty
DELTA = 30

# Mismatch cost matrix α(x, y)
ALPHA = {
    ('A', 'C'): 110,  ('A', 'G'):  48,  ('A', 'T'):  94,
    ('C', 'G'): 118,  ('C', 'T'):  48,
    ('G', 'T'): 110,
    # symmetric: α(x, y) == α(y, x), α(x, x) == 0
}
```

---

## Output

Each run prints:

| Field | Description |
|---|---|
| **Alignment cost** | Total optimal alignment cost |
| **Aligned X** | Sequence X with gap characters inserted |
| **Aligned Y** | Sequence Y with gap characters inserted |
| **Time (seconds)** | Wall-clock execution time |
| **Memory (MB)** | Peak process memory usage |

---

## Requirements

- Python 3.8+
- No external dependencies (standard library only)

---

## Background

Sequence alignment is a foundational problem in bioinformatics, used to measure similarity between DNA, RNA, or protein sequences. The Needleman–Wunsch algorithm (1970) introduced global alignment via dynamic programming. Hirschberg's algorithm (1975) later showed it could be solved in linear space — critical when aligning long genomic sequences.
