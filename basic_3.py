import sys
import time
import tracemalloc
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


def _basic_sequence_alignment(string1, string2):

    alpha = {'AA': 0, 'CC': 0, 'GG': 0, 'TT': 0, 'AC': 110, 'CA': 110, 'AG': 48, 'GA': 48, 'AT': 94, 'TA': 94,
            'CG': 118, 'GC': 118, 'CT': 48, 'TC': 48, 'GT': 110, 'TG': 110}
    delta = 30

    m,n = len(string1), len(string2)
    #initialize the dp array
    dp = [[0] * (n + 1) for i in range(m + 1)]

    #initialize the first row and column with gap penalties
    for i in range(0, n + 1):
        dp[0][i] = i * delta

    for i in range(1, m + 1):
        dp[i][0] = i * delta

    #fill the dp array
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            string = string1[i - 1] + string2[j - 1]
            dp[i][j] = min(
                dp[i - 1][j - 1] + alpha[string],
                dp[i-1][j] + delta, 
                dp[i][j-1] + delta
            )

    #backtrack to find the optimal alignment       
    i, j = m, n
    x_chars = []
    y_chars = []

    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + alpha[string1[i - 1] + string2[j - 1]]:
            x_chars.append(string1[i - 1])
            y_chars.append(string2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + delta:
            x_chars.append(string1[i - 1])
            y_chars.append("_")
            i -= 1
        else:
            x_chars.append("_")
            y_chars.append(string2[j - 1])
            j -= 1

    x = ''.join(reversed(x_chars))
    y = ''.join(reversed(y_chars))

    
    return dp[m][n], x, y


def main(argv: List[str]) -> None:
    if len(argv) != 3:
        print("usage: python csci570.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    in_path = Path(argv[1])
    out_path = Path(argv[2])

    s1, idx1, s2, idx2 = _parse_input(in_path)

    expanded1 = _expand(s1, idx1)
    expanded2 = _expand(s2, idx2)

    if len(expanded1) != (2 ** len(idx1)) * len(s1):
        raise RuntimeError("Length check failed for first string – verify indices.")
    if len(expanded2) != (2 ** len(idx2)) * len(s2):
        raise RuntimeError("Length check failed for second string – verify indices.")

    with out_path.open("w") as fh:
        fh.write(expanded1 + "\n")
        fh.write(expanded2 + "\n")

    #print("First expanded string :", expanded1)
    #print("Second expanded string:", expanded2)

    start =time.time()
    tracemalloc.start()
    Output = _basic_sequence_alignment(expanded1, expanded2)
    memory_usage = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end=time.time()
    runtime = (end-start)*1000

    f = open(out_path, 'w')
    f.write(str(int(Output[0])) + "\n")
    f.write(Output[1] + "\n")
    f.write(Output[2] + "\n")
    f.write(str(runtime) + "\n")
    f.write(str(memory_usage[1] / 1024) + "\n")
    f.close()


if __name__ == "__main__":
    main(sys.argv)
