"""Generate randomized, strength-2 Taguchi-style DOE tables.

Core engine used by the Streamlit app (app.py). It can also be used on its own:

    from taguchi_doe import generate_doe, print_report

    factors = {
        "Temperature [K]": [353.15, 358.15, 363.15, 368.15, 373.15],
        "Velocity [m/s]": [0.5, 1.0, 1.5, 2.0, 2.5],
        "Material": ["A", "B", "C", "D", "E"],
    }
    runs, report = generate_doe(factors, seed=42)
    print_report(report)

Matrix and degrees-of-freedom notes
-----------------------------------
For a q-level factor, the main effect uses q-1 degrees of freedom. An
orthogonal array OA(N, k, q, 2) has N rows and k q-level columns; strength 2
means every pair of columns contains every level pair equally often
(lambda = N/q**2). The prime-power-level arrays below are constructed over
finite fields with q symbols. The 50-run five-level array is the established
OA(50, 11, 5, 2) used for 7-11 five-level factors.

For a q-level group, the generated family uses N=q**m runs and up to
(q**m-1)/(q-1) columns; it selects the smallest m that can hold the factors.
This is an orthogonal-array capacity relation, not a guarantee that all
interactions are estimable.

The report shows N-1 total degrees of freedom, main-effect degrees of
freedom, and the remainder after fitting only main effects. That remainder
is not automatically pure experimental error: it can contain interactions
or other unmodeled effects. This generator does not assign interaction
columns. A two-factor interaction between factors with a and b levels needs
(a-1)*(b-1) degrees of freedom and should be checked for aliasing before
using the design to estimate it.

Taguchi mode supports factors with a prime-power number of levels (2, 3, 4,
5, 7, 8, 9, 11, 16, ...). Factors may have different supported level counts; their OA groups
are combined as a Cartesian product, preserving pairwise balance. For
non-prime-power level counts such as 6, 10, or 12, use
mode="full_factorial". The run
limit is checked before generating the table.
"""

from __future__ import annotations

import csv
import html
import itertools
import math
import random
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence


# OA(50, 11, 5, 2). The source array has one two-level column followed by
# eleven five-level columns; this generator uses only the five-level columns.
_L50_TEXT = """
1 1 1 1 1 1 1 1 1 1 1 1
1 1 2 2 2 2 2 2 2 2 2 2
1 1 3 3 3 3 3 3 3 3 3 3
1 1 4 4 4 4 4 4 4 4 4 4
1 1 5 5 5 5 5 5 5 5 5 5
1 2 1 2 3 4 5 1 2 3 4 5
1 2 2 3 4 5 1 2 3 4 5 1
1 2 3 4 5 1 2 3 4 5 1 2
1 2 4 5 1 2 3 4 5 1 2 3
1 2 5 1 2 3 4 5 1 2 3 4
1 3 1 3 5 2 4 4 1 3 5 2
1 3 2 4 1 3 5 5 2 4 1 3
1 3 3 5 2 4 1 1 3 5 2 4
1 3 4 1 3 5 2 2 4 1 3 5
1 3 5 2 4 1 3 3 5 2 4 1
1 4 1 4 2 5 3 5 3 1 4 2
1 4 2 5 3 1 4 1 4 2 5 3
1 4 3 1 4 2 5 2 5 3 1 4
1 4 4 2 5 3 1 3 1 4 2 5
1 4 5 3 1 4 2 4 2 5 3 1
1 5 1 5 4 3 2 4 3 2 1 5
1 5 2 1 5 4 3 5 4 3 2 1
1 5 3 2 1 5 4 1 5 4 3 2
1 5 4 3 2 1 5 2 1 5 4 3
1 5 5 4 3 2 1 3 2 1 5 4
2 1 1 1 4 5 4 3 2 5 2 3
2 1 2 2 5 1 5 4 3 1 3 4
2 1 3 3 1 2 1 5 4 2 4 5
2 1 4 4 2 3 2 1 5 3 5 1
2 1 5 5 3 4 3 2 1 4 1 2
2 2 1 2 1 3 3 2 4 5 5 4
2 2 2 3 2 4 4 3 5 1 1 5
2 2 3 4 3 5 5 4 1 2 2 1
2 2 4 5 4 1 1 5 2 3 3 2
2 2 5 1 5 2 2 1 3 4 4 3
2 3 1 3 3 1 2 5 5 4 2 4
2 3 2 4 4 2 3 1 1 5 3 5
2 3 3 5 5 3 4 2 2 1 4 1
2 3 4 1 1 4 5 3 3 2 5 2
2 3 5 2 2 5 1 4 4 3 1 3
2 4 1 4 5 4 1 2 5 2 3 3
2 4 2 5 1 5 2 3 1 3 4 4
2 4 3 1 2 1 3 4 2 4 5 5
2 4 4 2 3 2 4 5 3 5 1 1
2 4 5 3 4 3 5 1 4 1 2 2
2 5 1 5 2 2 5 3 4 4 3 1
2 5 2 1 3 3 1 4 5 5 4 2
2 5 3 2 4 4 2 5 1 1 5 3
2 5 4 3 5 5 3 1 2 2 1 4
2 5 5 4 1 1 4 2 3 3 2 5
"""
_L50 = [tuple(int(v) for v in line.split()) for line in _L50_TEXT.strip().splitlines()]

# Published L18 mixed array: one two-level column and seven three-level columns.
_L18_TEXT = """
1 1 1 1 1 1 1 1
1 1 2 2 2 2 2 2
1 1 3 3 3 3 3 3
1 2 1 1 2 2 3 3
1 2 2 2 3 3 1 1
1 2 3 3 1 1 2 2
1 3 1 2 1 3 2 3
1 3 2 3 2 1 3 1
1 3 3 1 3 2 1 2
2 1 1 3 3 2 2 1
2 1 2 1 1 3 3 2
2 1 3 2 2 1 1 3
2 2 1 2 3 1 3 2
2 2 2 3 1 2 1 3
2 2 3 1 2 3 2 1
2 3 1 3 2 3 1 2
2 3 2 1 3 1 2 3
2 3 3 2 1 2 3 1
"""
_L18 = [tuple(int(v) for v in line.split()) for line in _L18_TEXT.strip().splitlines()]


def _prime_power_parameters(q: int) -> tuple[int, int] | None:
    """Return (prime, exponent) when q is a prime power."""
    if q < 2:
        return None
    p = 2
    while p * p <= q:
        if q % p == 0:
            remainder = q
            exponent = 0
            while remainder % p == 0:
                remainder //= p
                exponent += 1
            return (p, exponent) if remainder == 1 else None
        p = 3 if p == 2 else p + 2
    return q, 1


def _poly_remainder(dividend: Sequence[int], divisor: Sequence[int], p: int) -> tuple[int, ...]:
    """Polynomial remainder over GF(p), with coefficients in ascending order."""
    remainder = [coefficient % p for coefficient in dividend]
    while len(remainder) > 1 and remainder[-1] == 0:
        remainder.pop()
    while len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        leading = remainder[-1]  # each divisor passed here is monic
        for i, coefficient in enumerate(divisor):
            remainder[shift + i] = (remainder[shift + i] - leading * coefficient) % p
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return tuple(remainder)


def _irreducible_polynomial(p: int, degree: int) -> tuple[int, ...]:
    """Find a monic irreducible polynomial of the requested degree over GF(p)."""
    if degree == 1:
        return (0, 1)
    for lower_coefficients in itertools.product(range(p), repeat=degree):
        if lower_coefficients[0] == 0:
            continue
        candidate = tuple(lower_coefficients) + (1,)
        reducible = False
        for factor_degree in range(1, degree // 2 + 1):
            for lower in itertools.product(range(p), repeat=factor_degree):
                divisor = tuple(lower) + (1,)
                if _poly_remainder(candidate, divisor, p) == (0,):
                    reducible = True
                    break
            if reducible:
                break
        if not reducible:
            return candidate
    raise ValueError(f"Could not construct GF({p}^{degree}) for {p**degree} levels.")


def _field_tables(q: int) -> tuple[list[list[int]], list[list[int]]]:
    """Create addition and multiplication tables for a finite field of size q."""
    parameters = _prime_power_parameters(q)
    if parameters is None:
        raise ValueError(f"{q} is not a prime-power field size.")
    p, degree = parameters
    modulus = _irreducible_polynomial(p, degree)

    def digits(value: int) -> list[int]:
        result = []
        for _ in range(degree):
            result.append(value % p)
            value //= p
        return result

    def encode(coefficients: Sequence[int]) -> int:
        return sum((coefficient % p) * (p**i) for i, coefficient in enumerate(coefficients))

    additions = [[0] * q for _ in range(q)]
    multiplications = [[0] * q for _ in range(q)]
    digit_cache = [digits(value) for value in range(q)]
    for left in range(q):
        a = digit_cache[left]
        for right in range(q):
            b = digit_cache[right]
            additions[left][right] = encode([(x + y) % p for x, y in zip(a, b)])
            product_coefficients = [0] * (2 * degree - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    product_coefficients[i + j] += x * y
            product_coefficients = [x % p for x in product_coefficients]
            for power in range(len(product_coefficients) - 1, degree - 1, -1):
                leading = product_coefficients[power] % p
                if leading:
                    shift = power - degree
                    for i in range(degree):
                        product_coefficients[shift + i] = (
                            product_coefficients[shift + i] - leading * modulus[i]
                        ) % p
            multiplications[left][right] = encode(product_coefficients[:degree])
    return additions, multiplications


def _validate_oa(rows: Sequence[Sequence[int]], q: int, columns: int) -> None:
    if not rows or any(len(row) != columns for row in rows):
        raise ValueError("Internal DOE array has an invalid shape.")
    n_runs = len(rows)
    for col in range(columns):
        counts = Counter(row[col] for row in rows)
        if set(counts) != set(range(q)) or set(counts.values()) != {n_runs // q}:
            raise ValueError(f"Internal DOE column {col + 1} is not balanced.")
    if columns > 1:
        if n_runs % (q * q):
            raise ValueError("Internal array run count cannot support strength 2.")
        expected = n_runs // (q * q)
        for left in range(columns):
            for right in range(left + 1, columns):
                counts = Counter((row[left], row[right]) for row in rows)
                if len(counts) != q * q or set(counts.values()) != {expected}:
                    raise ValueError(
                        f"Internal DOE columns {left + 1} and {right + 1} "
                        "are not pairwise orthogonal."
                    )


def _make_oa(
    q: int, n_factors: int, max_runs: int, max_cells: int, extra_columns: int
) -> tuple[list[tuple[int, ...]], dict[str, Any]]:
    """Build a strength-2 OA for a prime-power q using projective vectors."""
    if q == 5 and 7 <= n_factors <= 11:
        if len(_L50) * (2 * n_factors + extra_columns) > max_cells:
            raise ValueError(f"The selected five-level OA exceeds max_cells={max_cells:,}.")
        rows = [tuple(value - 1 for value in row[1 : n_factors + 1]) for row in _L50]
        _validate_oa(rows, 5, n_factors)
        return rows, {
            "array": "L50 / OA(50, 11, 5, 2)",
            "columns_used": list(range(2, n_factors + 2)),
            "runs": len(rows),
            "pair_frequency": 2,
        }

    dimension = 1
    while (q**dimension - 1) // (q - 1) < n_factors:
        dimension += 1
        if q**dimension > max_runs:
            max_columns = (q**dimension - 1) // (q - 1)
            raise ValueError(
                f"Need {n_factors} factors at {q} levels. The next supported "
                f"prime-power OA needs {q**dimension:,} runs for up to "
                f"{max_columns} columns, above max_runs={max_runs:,}."
            )

    n_runs = q**dimension
    if n_runs > max_runs:
        raise ValueError(
            f"The selected {q}-level OA needs {n_runs:,} runs, above "
            f"max_runs={max_runs:,}."
        )
    if n_runs * (2 * n_factors + extra_columns) > max_cells:
        raise ValueError(
            f"The selected {q}-level OA exceeds max_cells={max_cells:,} "
            "for its run and factor count."
        )

    if n_factors == 1:
        rows = [(level,) for level in range(q)]
        _validate_oa(rows, q, 1)
        return rows, {
            "array": f"L{q} single-factor design",
            "columns_used": [1],
            "runs": n_runs,
            "pair_frequency": None,
        }

    additions, multiplications = _field_tables(q)
    # One representative from each nonzero scalar-multiple class gives the
    # (q**dimension - 1)/(q - 1) columns of the projective OA.
    column_vectors = []
    for vector in itertools.product(range(q), repeat=dimension):
        first_nonzero = next((value for value in vector if value != 0), None)
        if first_nonzero == 1:
            column_vectors.append(vector)
    all_rows = itertools.product(range(q), repeat=dimension)
    rows = []
    for row_vector in all_rows:
        row = []
        for col_vector in column_vectors[:n_factors]:
            value = 0
            for a, b in zip(row_vector, col_vector):
                value = additions[value][multiplications[a][b]]
            row.append(value)
        rows.append(tuple(row))
    _validate_oa(rows, q, n_factors)
    array_name = f"OA({n_runs}, {(q**dimension - 1) // (q - 1)}, {q}, 2)"
    return rows, {
        "array": array_name,
        "columns_used": list(range(1, n_factors + 1)),
        "runs": n_runs,
        "pair_frequency": n_runs // (q * q) if n_factors > 1 else None,
    }


def _validate_mixed_oa(
    rows: Sequence[Mapping[str, int]], factor_names: Sequence[str], level_counts: Mapping[str, int]
) -> None:
    """Check equal marginal and pair frequencies in a coded mixed-level array."""
    n_runs = len(rows)
    for name in factor_names:
        q = level_counts[name]
        if n_runs % q:
            raise ValueError(f"Internal design is not balanced for factor {name!r}.")
        counts = Counter(row[name] for row in rows)
        if set(counts) != set(range(q)) or set(counts.values()) != {n_runs // q}:
            raise ValueError(f"Internal design is not balanced for factor {name!r}.")
    for i, left in enumerate(factor_names):
        for right in factor_names[i + 1 :]:
            pair_count = level_counts[left] * level_counts[right]
            if n_runs % pair_count:
                raise ValueError(f"Internal design cannot balance {left!r} and {right!r}.")
            expected = n_runs // pair_count
            counts = Counter((row[left], row[right]) for row in rows)
            if len(counts) != pair_count or set(counts.values()) != {expected}:
                raise ValueError(f"Internal design is not orthogonal for {left!r}, {right!r}.")


def _factorial_codes(factor_names: Sequence[str], level_counts: Mapping[str, int]) -> list[dict[str, int]]:
    combinations = itertools.product(*(range(level_counts[name]) for name in factor_names))
    return [dict(zip(factor_names, combination)) for combination in combinations]


def generate_doe(
    factors: Mapping[str, Sequence[Any]],
    *,
    mode: str = "taguchi",
    response_columns: Sequence[str] = (),
    include_codes: bool = False,
    randomize: bool = True,
    seed: int | None = 42,
    replicates: int = 1,
    max_runs: int = 100_000,
    max_cells: int = 2_000_000,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Generate DOE run records and a design/DOF report.

    Parameters
    ----------
    factors:
        Ordered mapping of factor names to their actual numeric or categorical
        levels. Level order is preserved and codes start at 1.
    mode:
        ``"taguchi"`` uses supported prime-power-level strength-2 arrays.
        ``"full_factorial"`` enumerates all combinations and accepts any number
        of levels.
    response_columns:
        Optional names for blank result columns to fill after numerical, CFD,
        or physical runs.
    include_codes:
        Include coded level columns in the returned records. The default table
        view shows actual factor settings only, like a lab run sheet.
    randomize:
        Randomize physical execution order after constructing standard order.
    seed:
        Reproducible randomization seed; use ``None`` for a new order each run.
    replicates:
        Number of repeats of every design point. Repeats are separately
        randomized when ``randomize`` is true.
    max_runs:
        Maximum output rows, including replicates, to prevent accidental
        generation of very large designs.
    max_cells:
        Maximum approximate table cells (rows x columns), including factor
        columns and any optional coded levels.

    Each record contains run order, standard run, replicate, actual factor
    values, optional 1-based coded levels (``code:<factor>``), and any blank
    response columns requested.
    """
    if not isinstance(factors, Mapping) or not factors:
        raise ValueError("factors must be a non-empty mapping of names to levels.")
    if mode not in {"taguchi", "full_factorial"}:
        raise ValueError("mode must be 'taguchi' or 'full_factorial'.")
    if not isinstance(randomize, bool):
        raise ValueError("randomize must be True or False.")
    if not isinstance(include_codes, bool):
        raise ValueError("include_codes must be True or False.")
    if not isinstance(replicates, int) or isinstance(replicates, bool) or replicates < 1:
        raise ValueError("replicates must be a positive integer.")
    if not isinstance(max_runs, int) or max_runs < 1:
        raise ValueError("max_runs must be a positive integer.")
    if not isinstance(max_cells, int) or max_cells < 1:
        raise ValueError("max_cells must be a positive integer.")

    factor_names = list(factors.keys())
    reserved = {"Run", "RunOrder", "StandardRun", "Replicate", "Status"}
    if any(not isinstance(name, str) or not name.strip() for name in factor_names):
        raise ValueError("Every factor name must be a non-empty string.")
    if any(name in reserved or name.startswith("code:") for name in factor_names):
        raise ValueError("Factor names cannot use RunOrder, StandardRun, Replicate, or the 'code:' prefix.")
    response_names = list(response_columns)
    if any(not isinstance(name, str) or not name.strip() for name in response_names):
        raise ValueError("Every response column name must be a non-empty string.")
    if len(set(response_names)) != len(response_names):
        raise ValueError("Response column names must be unique.")
    if set(response_names) & (set(factor_names) | reserved) or any(name.startswith("code:") for name in response_names):
        raise ValueError("Response columns cannot duplicate factor, reserved, or coded-level column names.")

    levels: dict[str, list[Any]] = {}
    level_counts: dict[str, int] = {}
    for name in factor_names:
        values = list(factors[name])
        if len(values) < 2:
            raise ValueError(f"Factor {name!r} must have at least two levels.")
        for value in values:
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"Factor {name!r} contains a non-finite numeric level.")
        try:
            duplicate_levels = len(set(values)) != len(values)
        except TypeError:
            duplicate_levels = any(
                values[i] == values[j] for i in range(len(values)) for j in range(i)
            )
        if duplicate_levels:
            raise ValueError(f"Factor {name!r} contains duplicate levels.")
        levels[name] = values
        level_counts[name] = len(values)

    design_groups: list[dict[str, Any]] = []
    output_columns = 3 + len(factor_names) * (2 if include_codes else 1) + len(response_names)
    if mode == "full_factorial":
        base_n = math.prod(level_counts.values())
        if base_n * replicates > max_runs:
            raise ValueError(
                f"Full factorial needs {base_n * replicates:,} rows including "
                f"replicates, above max_runs={max_runs:,}."
            )
        if base_n * replicates * output_columns > max_cells:
            raise ValueError(f"Full factorial table exceeds max_cells={max_cells:,}.")
        coded_rows = _factorial_codes(factor_names, level_counts)
        method = "full factorial"
    else:
        groups: dict[int, list[str]] = {}
        for name in factor_names:
            q = level_counts[name]
            if _prime_power_parameters(q) is None:
                raise ValueError(
                    f"Factor {name!r} has {q} levels. Taguchi mode currently "
                    "supports prime-power level counts; use mode='full_factorial' "
                    "for this factor or change to a supported level count."
                )
            groups.setdefault(q, []).append(name)

        if set(groups) == {2, 3} and len(groups[2]) == 1 and 4 <= len(groups[3]) <= 7:
            # The L18 is smaller than the Cartesian product for one binary and
            # four-to-seven ternary factors. Validate the selected columns below.
            binary_name = groups[2][0]
            ternary_names = groups[3]
            coded_rows = [
                {binary_name: row[0] - 1,
                 **{name: row[i + 1] - 1 for i, name in enumerate(ternary_names)}}
                for row in _L18
            ]
            base_n = len(coded_rows)
            design_groups.append({
                "level_count": "mixed 2/3",
                "factors": [binary_name, *ternary_names],
                "array": "L18 / OA(18, 1x2 + 7x3, 2)",
                "columns_used": [1, *range(2, len(ternary_names) + 2)],
                "runs": base_n,
                "pair_frequency": "3 for 2x3; 2 for 3x3",
            })
            method = "mixed-level L18 Taguchi orthogonal array"
        else:
            group_matrices: list[tuple[list[str], list[tuple[int, ...]]]] = []
            for q, names in groups.items():
                matrix, info = _make_oa(q, len(names), max_runs, max_cells, 3 + len(response_names))
                design_groups.append({"level_count": q, "factors": names, **info})
                group_matrices.append((names, matrix))

            base_n = math.prod(len(matrix) for _, matrix in group_matrices)
            if base_n * replicates > max_runs:
                raise ValueError(
                    f"The combined mixed-level Taguchi design needs "
                    f"{base_n * replicates:,} rows including replicates, above "
                    f"max_runs={max_runs:,}."
                )
            if base_n * replicates * output_columns > max_cells:
                raise ValueError(
                    f"Combined Taguchi table exceeds max_cells={max_cells:,}."
                )
            coded_rows = []
            for group_combination in itertools.product(*(matrix for _, matrix in group_matrices)):
                coded: dict[str, int] = {}
                for (names, _), group_row in zip(group_matrices, group_combination):
                    coded.update(dict(zip(names, group_row)))
                coded_rows.append(coded)
            method = "prime-power strength-2 Taguchi orthogonal array"

        if base_n * replicates > max_runs:
            raise ValueError(
                f"The selected Taguchi design needs {base_n * replicates:,} rows "
                f"including replicates, above max_runs={max_runs:,}."
            )

    if len(coded_rows) != base_n:
        raise ValueError("Internal DOE row count does not match the selected design.")
    _validate_mixed_oa(coded_rows, factor_names, level_counts)

    main_effect_dof = sum(q - 1 for q in level_counts.values())
    base_total_dof = base_n - 1
    if main_effect_dof > base_total_dof:
        raise ValueError(
            f"Main effects require {main_effect_dof} DOF, but this design has "
            f"only {base_total_dof} total DOF."
        )

    records: list[dict[str, Any]] = []
    for replicate in range(1, replicates + 1):
        for standard_run, coded in enumerate(coded_rows, start=1):
            record: dict[str, Any] = {
                "RunOrder": None,
                "StandardRun": standard_run,
                "Replicate": replicate,
            }
            for name in factor_names:
                if include_codes:
                    record[f"code:{name}"] = coded[name] + 1
                record[name] = levels[name][coded[name]]
            for response_name in response_names:
                record[response_name] = ""
            records.append(record)

    if randomize:
        random.Random(seed).shuffle(records)
    for run_order, record in enumerate(records, start=1):
        record["RunOrder"] = run_order

    total_runs = len(records)
    total_dof = total_runs - 1
    total_main_effect_dof = main_effect_dof
    report = {
        "method": method,
        "factors": len(factor_names),
        "factor_level_counts": level_counts,
        "response_columns": response_names,
        "design_groups": design_groups,
        "base_design_runs": base_n,
        "replicates": replicates,
        "total_runs": total_runs,
        "randomized": randomize,
        "random_seed": seed if randomize else None,
        "strength": 2 if len(factor_names) > 1 else 1,
        "main_effect_dof": total_main_effect_dof,
        "total_dof_after_replicates": total_dof,
        "remaining_dof_after_main_effects": total_dof - total_main_effect_dof,
        "base_design_total_dof": base_total_dof,
        "base_design_remaining_dof_after_main_effects": base_total_dof - main_effect_dof,
        "factor_dof": {name: level_counts[name] - 1 for name in factor_names},
    }
    return records, report


def print_report(report: Mapping[str, Any]) -> None:
    """Print the design selection and main-effect degrees-of-freedom summary."""
    print(f"Design: {report['method']}")
    print(f"Factors: {report['factors']}; runs: {report['total_runs']} "
          f"({report['base_design_runs']} unique design points x {report['replicates']} replicate(s))")
    print(f"Strength: {report['strength']}")
    for group in report["design_groups"]:
        print(f"  {group['array']}: {len(group['factors'])} factor(s), "
              f"columns {group['columns_used']}, pair frequency {group['pair_frequency']}")
    factor_dof = ", ".join(f"{name}: {dof}" for name, dof in report["factor_dof"].items())
    print(f"Main-effect DOF: {report['main_effect_dof']} ({factor_dof})")
    print(f"Base design total DOF (N-1): {report['base_design_total_dof']}")
    print("Base design DOF left after main effects: "
          f"{report['base_design_remaining_dof_after_main_effects']}")
    print(f"Total DOF after replication (N-1): {report['total_dof_after_replicates']}")
    print("DOF left after main effects including repeats: "
          f"{report['remaining_dof_after_main_effects']}")
    if report["randomized"]:
        print(f"Run order randomized with seed {report['random_seed']}.")
    print("Remaining DOF can include interactions or other unmodeled effects; it is not automatically pure error.")


def evaluate_doe(
    records: Sequence[Mapping[str, Any]],
    factor_names: Sequence[str],
    response_columns: Sequence[str],
    evaluator: Any,
) -> list[dict[str, Any]]:
    """Run a user-supplied model over every DOE run and return a results table.

    ``evaluator`` receives a mapping of actual factor settings and must return
    a mapping from each response column name to its result. Exceptions are
    recorded in ``Status`` so one invalid run does not hide the rest.
    """
    if not records:
        raise ValueError("records is empty; there are no DOE runs to evaluate.")
    if not callable(evaluator):
        raise ValueError("evaluator must be a callable model function.")
    results = []
    for record in records:
        case = {name: record[name] for name in factor_names}
        row: dict[str, Any] = {"Run": record["RunOrder"], **case}
        try:
            output = evaluator(case)
            if not isinstance(output, Mapping):
                raise TypeError("model evaluator must return a mapping of response names to values")
            missing = [name for name in response_columns if name not in output]
            if missing:
                raise ValueError(f"model did not return response(s): {', '.join(missing)}")
            row.update({name: output[name] for name in response_columns})
            row["Status"] = "OK"
        except Exception as exc:
            row.update({name: "" for name in response_columns})
            row["Status"] = f"INVALID: {exc}"
        results.append(row)
    return results


def write_csv(records: Sequence[Mapping[str, Any]], path: str | Path) -> None:
    """Write a generated DOE table as UTF-8 CSV, with Excel-friendly BOM."""
    if not records:
        raise ValueError("records is empty; there is no DOE table to write.")
    fieldnames = list(records[0].keys())
    with Path(path).open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
