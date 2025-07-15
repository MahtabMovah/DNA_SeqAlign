import time
import tracemalloc
import sys
import sys
from pathlib import Path
from typing import List, Tuple

def _parse_input(path: Path) -> Tuple[str, List[int], str, List[int]]:

    with path.open() as fh:
        raw = [ln.strip() for ln in fh if ln.strip()]

    if not raw:
        raise ValueError("Input file is empty.")

    # 1) first base string
    s1 = raw[0]

    # 2) indices for the first string – stop at the first non‑numeric token
    idx1: List[int] = []
    i = 1
    while i < len(raw) and raw[i].isdigit():
        idx1.append(int(raw[i]))
        i += 1

    if i >= len(raw):
        raise ValueError("Second base string is missing.")

    # 3) second base string
    s2 = raw[i]
    i += 1

    # 4) indices for the second string – the rest of the file must be numeric
    idx2: List[int] = []
    while i < len(raw):
        if not raw[i].isdigit():
            raise ValueError(f"Unexpected non‑numeric token: {raw[i]}")
        idx2.append(int(raw[i]))
        i += 1

    return s1, idx1, s2, idx2


def _expand(base: str, indices: List[int]) -> str:
    """Iteratively expand *base* using *indices* (0‑based, inclusive)."""
    s = base
    for idx in indices:
        if idx < 0 or idx >= len(s):
            raise IndexError(
                f"Index {idx} is out of bounds for current string length {len(s)}"
            )
        prev = s  # copy of the current cumulative string
        # insert the copy *after* character at position idx
        s = prev[: idx + 1] + prev + prev[idx + 1 :]
    return s

# ---- mismatch cost matrix and gap penalty ----
δ_e = 30
α = {
    'A': {'A': 0,   'C': 110, 'G': 48,  'T': 94},
    'C': {'A': 110, 'C': 0,   'G': 118, 'T': 48},
    'G': {'A': 48,  'C': 118, 'G': 0,   'T': 110},
    'T': {'A': 94,  'C': 48,  'G': 110, 'T': 0}
}

# ---- helper: Needleman-Wunsch for small cases ----
def needleman_wunsch(X, Y, gap, alpha):
    m, n = len(X), len(Y)
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(1, m+1):
        dp[i][0] = i * gap
    for j in range(1, n+1):
        dp[0][j] = j * gap

    for i in range(1, m+1):
        xi = X[i-1]
        for j in range(1, n+1):
            yj = Y[j-1]
            dp[i][j] = min(
                dp[i-1][j-1] + alpha[xi][yj],
                dp[i-1][j]   + gap,
                dp[i][j-1]   + gap
            )

    # backtrack
    A, B = [], []
    i, j = m, n
    while i>0 or j>0:
        if i>0 and j>0 and dp[i][j] == dp[i-1][j-1] + alpha[X[i-1]][Y[j-1]]:
            A.append(X[i-1]);  B.append(Y[j-1])
            i -= 1; j -= 1
        elif i>0 and dp[i][j] == dp[i-1][j] + gap:
            A.append(X[i-1]);  B.append('-')
            i -= 1
        else:
            A.append('-');     B.append(Y[j-1])
            j -= 1

    return ''.join(reversed(A)), ''.join(reversed(B))


def nw_score(X, Y, gap, alpha):
    n = len(Y)
    prev = [j*gap for j in range(n+1)]
    for xi in X:
        curr = [prev[0] + gap] + [0]*n
        for j, yj in enumerate(Y, start=1):
            curr[j] = min(
                prev[j-1] + alpha[xi][yj],
                prev[j]   + gap,
                curr[j-1] + gap
            )
        prev = curr
    return prev


def hirschberg(X, Y, gap, alpha):
    m, n = len(X), len(Y)
    if m == 0:
        return '-'*n, Y
    if n == 0:
        return X, '-'*m
    if m == 1 or n == 1:
        return needleman_wunsch(X, Y, gap, alpha)

    i_mid = m // 2
    scoreL = nw_score(X[:i_mid], Y,       gap, alpha)
    scoreR = nw_score(X[i_mid:][::-1], Y[::-1], gap, alpha)
    # find split of Y that minimizes total cost
    total = [l + r for l, r in zip(scoreL, reversed(scoreR))]
    j_mid = total.index(min(total))

    A1, B1 = hirschberg(X[:i_mid], Y[:j_mid], gap, alpha)
    A2, B2 = hirschberg(X[i_mid:], Y[j_mid:], gap, alpha)
    return A1 + A2, B1 + B2

def align_sequences(X, Y, gap, alpha):
    tracemalloc.start()
    t0 = time.perf_counter()
    A, B = hirschberg(X, Y, gap, alpha)
    t1 = time.perf_counter()

    # compute final cost
    cost = 0
    for a, b in zip(A, B):
        if a == '-' or b == '-':
            cost += gap
        else:
            cost += alpha[a][b]

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    runtime_ms = (t1 - t0) * 1000.0
    peak_kb    = peak / 1024.0

    return cost, A, B, runtime_ms, peak_kb



def main(argv: List[str]) -> None:
    if len(argv) != 3:
        print("usage: python csci570.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    in_path = Path(argv[1])
    out_path = Path(argv[2])

    s1, idx1, s2, idx2 = _parse_input(in_path)

    expanded1 = _expand(s1, idx1)
    expanded2 = _expand(s2, idx2)

    #print("First expanded string :", expanded1)
    #print("Second expanded string:", expanded2)
    X =expanded1
    #print("Enter second sequence:")
    Y =expanded2
    
    cost, A, B, runtime_ms, peak_kb = align_sequences(X, Y, δ_e, α)

    if len(expanded1) != (2 ** len(idx1)) * len(s1):
        raise RuntimeError("Length check failed for first string – verify indices.")
    if len(expanded2) != (2 ** len(idx2)) * len(s2):
        raise RuntimeError("Length check failed for second string – verify indices.")

    with out_path.open("w") as fh:
        fh.write(str(cost) + "\n")
        fh.write(A + "\n")
        fh.write(B + "\n")
        fh.write(str(runtime_ms) + "\n")
        fh.write(str(peak_kb) + "\n")
        





    #print("\n=== Alignment results ===")
    #print(f"Total cost       : {cost}")
    #print(f"Aligned sequence1: {A}")
    #print(f"Aligned sequence2: {B}")
    #print(f"Runtime (ms)     : {runtime_ms:.2f}")
    #print(f"Peak memory (KB) : {peak_kb:.1f}")


if __name__ == "__main__":
    main(sys.argv)
