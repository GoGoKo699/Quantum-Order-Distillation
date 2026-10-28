"""Exact matrix-unit certificates for the qubit central-degree lower bound.

Multilinearity reduces scalar-valuedness to basis tuples.  For degrees 1,2,3,
a square subsystem of those linear conditions already has nonzero determinant.
The resulting certificate rules out every complex coefficient vector, not only
sampled polynomials.  Python standard library only.
"""
from fractions import Fraction
from itertools import permutations, product
from math import factorial
import json
from pathlib import Path

BASIS = [
    ((1, 0), (0, 0)),
    ((0, 1), (0, 0)),
    ((0, 0), (1, 0)),
    ((0, 0), (0, 1)),
]
COORDINATES = ("M01", "M10", "M00-M11")

def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))

def evaluate_word(word, values):
    out = ((1, 0), (0, 1))
    for j in word:
        out = matmul(out, values[j])
    return out

def scalar_obstructions(a):
    return a[0][1], a[1][0], a[0][0] - a[1][1]

def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    det = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            det = -det
        value = a[j][j]
        det *= value
        a[j] = [x / value for x in a[j]]
        for i in range(j + 1, len(a)):
            factor = a[i][j]
            if factor:
                a[i] = [x - factor*y for x, y in zip(a[i], a[j])]
    return det

def certificate(degree):
    words = list(permutations(range(degree)))
    pivots = {}
    rows = []
    labels = []
    for indices in product(range(4), repeat=degree):
        values = [BASIS[i] for i in indices]
        obs = [scalar_obstructions(evaluate_word(word, values)) for word in words]
        for coord in range(3):
            original = [x[coord] for x in obs]
            row = [Fraction(x) for x in original]
            for col in sorted(pivots):
                factor = row[col]
                if factor:
                    row = [x - factor*y for x, y in zip(row, pivots[col])]
            pivot = next((j for j, x in enumerate(row) if x), None)
            if pivot is None:
                continue
            scale = row[pivot]
            pivots[pivot] = [x / scale for x in row]
            rows.append(original)
            labels.append({"matrix_unit_indices": list(indices),
                           "coordinate": COORDINATES[coord]})
            if len(rows) == len(words):
                det = determinant(rows)
                assert det != 0
                return {"degree": degree, "coefficient_count": factorial(degree),
                        "rank": len(rows), "determinant": str(det),
                        "word_order": [list(w) for w in words],
                        "row_labels": labels, "certificate_matrix": rows}
    raise AssertionError("Did not obtain full rank")

def main():
    data = {"convention": "E00,E01,E10,E11; words multiply left to right",
            "certificates": [certificate(m) for m in (1, 2, 3)]}
    output = Path(__file__).with_name("qubit_minimum_results.json")
    expected = json.loads(output.read_text()) if output.exists() else None
    if expected is not None:
        assert data == expected, "Exact certificate differs from committed result"
    else:
        output.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({"qubit_minimum_certificates": [
        {"degree": x["degree"], "rank": x["rank"],
         "determinant": x["determinant"]} for x in data["certificates"]]}))

if __name__ == "__main__":
    main()

