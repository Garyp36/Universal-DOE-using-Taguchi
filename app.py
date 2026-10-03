"""Streamlit front end for the Taguchi DOE generator."""

from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st

from taguchi_doe import _prime_power_parameters, generate_doe

st.set_page_config(page_title="Taguchi DOE Generator", page_icon="🧪", layout="wide")

DEFAULT_FACTORS = pd.DataFrame(
    {
        "Factor": ["Temperature [K]", "Velocity [m/s]", "Material"],
        "Levels (comma-separated)": [
            "353.15, 358.15, 363.15, 368.15, 373.15",
            "0.5, 1.0, 1.5, 2.0, 2.5",
            "A, B, C, D, E",
        ],
    }
)


# ---------------------------------------------------------------- helpers
def parse_levels(text: str) -> list[Any]:
    """Split a comma-separated string. All-numeric -> numbers, otherwise text."""
    tokens = [t.strip() for t in str(text).split(",") if t.strip()]
    try:
        values: list[Any] = []
        for token in tokens:
            number = float(token)
            is_plain_int = number.is_integer() and not any(c in token for c in ".eE")
            values.append(int(number) if is_plain_int else number)
        return values
    except ValueError:
        return tokens


def read_factors(table: pd.DataFrame) -> tuple[dict[str, list[Any]], list[str]]:
    """Convert the editor table to a factor dict; also return problems found."""
    factors: dict[str, list[Any]] = {}
    problems: list[str] = []
    for position, row in table.iterrows():
        name = str(row.get("Factor") or "").strip()
        raw = row.get("Levels (comma-separated)")
        levels_text = "" if raw is None or (isinstance(raw, float) and pd.isna(raw)) else str(raw)
        if not name and not levels_text.strip():
            continue  # fully blank row
        label = name or f"row {position + 1}"
        if not name:
            problems.append(f"Row {position + 1} has levels but no factor name.")
            continue
        if name in factors:
            problems.append(f"Factor name '{name}' appears more than once.")
            continue
        levels = parse_levels(levels_text)
        if len(levels) < 2:
            problems.append(f"'{label}' needs at least two levels.")
            continue
        if len(set(levels)) != len(levels):
            problems.append(f"'{label}' has duplicate levels.")
            continue
        factors[name] = levels
    return factors, problems


def factor_summary(factors: dict[str, list[Any]]) -> pd.DataFrame:
    rows = []
    for name, levels in factors.items():
        q = len(levels)
        rows.append(
            {
                "Factor": name,
                "Levels": q,
                "Main-effect DOF": q - 1,
                "Taguchi-compatible": "Yes" if _prime_power_parameters(q) else "No",
                "Values": ", ".join(str(v) for v in levels),
            }
        )
    return pd.DataFrame(rows)


def build_excel(runs: pd.DataFrame, summary: pd.DataFrame, factors: pd.DataFrame) -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        runs.to_excel(writer, sheet_name="Run sheet", index=False)
        factors.to_excel(writer, sheet_name="Factors", index=False)
        summary.to_excel(writer, sheet_name="Design summary", index=False)
    return buffer.getvalue()


def report_to_frame(report: dict[str, Any]) -> pd.DataFrame:
    items = [
        ("Method", report["method"]),
        ("Factors", report["factors"]),
        ("Unique design points", report["base_design_runs"]),
        ("Replicates", report["replicates"]),
        ("Total runs", report["total_runs"]),
        ("Strength", report["strength"]),
        ("Randomized", "Yes" if report["randomized"] else "No"),
        ("Random seed", report["random_seed"] if report["randomized"] else "n/a"),
        ("Main-effect DOF", report["main_effect_dof"]),
        ("Total DOF (N-1, with replicates)", report["total_dof_after_replicates"]),
        ("DOF left after main effects", report["remaining_dof_after_main_effects"]),
        ("Base design DOF (N-1)", report["base_design_total_dof"]),
        ("Base design DOF left after main effects", report["base_design_remaining_dof_after_main_effects"]),
    ]
    return pd.DataFrame(items, columns=["Item", "Value"]).astype(str)


# ---------------------------------------------------------------- header
st.title("🧪 Taguchi DOE Generator")
st.write(
    "Build a randomized, strength-2 orthogonal-array run sheet from your factors and levels. "
    "Add factors below, pick the options in the sidebar, and press **Generate**."
)

# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Design options")
    mode_label = st.radio(
        "Design type",
        ["Taguchi orthogonal array", "Full factorial"],
        help=(
            "Taguchi needs a prime-power number of levels per factor (2, 3, 4, 5, 7, 8, 9, 11, 16...). "
            "Use full factorial for 6, 10, 12 levels, at the cost of many more runs."
        ),
    )
    mode = "taguchi" if mode_label.startswith("Taguchi") else "full_factorial"

    replicates = st.number_input("Replicates per design point", min_value=1, max_value=50, value=1, step=1)
    randomize = st.checkbox("Randomize run order", value=True)
    seed_text = st.text_input(
        "Random seed",
        value="42",
        disabled=not randomize,
        help="Same seed gives the same order. Type 'none' for a fresh order each time.",
    )
    include_codes = st.checkbox("Include coded levels (1, 2, 3...)", value=False)

    st.subheader("Response columns")
    responses_text = st.text_input(
        "Blank result columns (optional)",
        placeholder="Q [W], Outlet T [K], dP [Pa]",
        help="Comma-separated. Added as empty columns to fill after each CFD or physical run.",
    )

    with st.expander("Safety limits"):
        max_runs = st.number_input("Maximum runs", min_value=1, value=5000, step=500)
        max_cells = st.number_input("Maximum table cells", min_value=1, value=500_000, step=50_000)

# ---------------------------------------------------------------- factor input
st.subheader("1. Define factors")
st.caption(
    "One row per factor. Type levels separated by commas, in the order you want them coded. "
    "Numbers and text both work. Use the + at the bottom of the table to add rows."
)

table = st.data_editor(
    DEFAULT_FACTORS,
    num_rows="dynamic",
    width="stretch",
    hide_index=True,
    key="factor_editor",
    column_config={
        "Factor": st.column_config.TextColumn("Factor", help="Include units, e.g. Velocity [m/s]", required=False),
        "Levels (comma-separated)": st.column_config.TextColumn(
            "Levels (comma-separated)", help="e.g. 0.5, 1.0, 1.5, 2.0, 2.5", required=False
        ),
    },
)

factors, problems = read_factors(table)
for message in problems:
    st.warning(message)

if factors:
    summary_df = factor_summary(factors)
    st.dataframe(summary_df, width="stretch", hide_index=True)
    incompatible = summary_df[summary_df["Taguchi-compatible"] == "No"]["Factor"].tolist()
    full_size = 1
    for levels in factors.values():
        full_size *= len(levels)
    if mode == "taguchi" and incompatible:
        st.error(
            "Not usable with Taguchi: " + ", ".join(f"**{n}**" for n in incompatible)
            + ". Change those to 2, 3, 4, 5, 7, 8, 9 or 11 levels, or switch to Full factorial in the sidebar."
        )
    st.caption(f"A full factorial of these factors would need {full_size:,} runs.")

# ---------------------------------------------------------------- generate
st.subheader("2. Generate")
if st.button("Generate DOE table", type="primary", disabled=not factors):
    try:
        seed: int | None
        if not randomize:
            seed = None
        elif seed_text.strip().lower() == "none":
            seed = None
        else:
            seed = int(seed_text.strip() or "42")
        response_columns = [p.strip() for p in responses_text.split(",") if p.strip()]
        runs, report = generate_doe(
            factors,
            mode=mode,
            response_columns=response_columns,
            include_codes=include_codes,
            randomize=randomize,
            seed=seed,
            replicates=int(replicates),
            max_runs=int(max_runs),
            max_cells=int(max_cells),
        )
        st.session_state["result"] = {
            "runs": runs,
            "report": report,
            "factors": factors,
            "response_columns": response_columns,
            "include_codes": include_codes,
        }
    except ValueError as exc:
        st.session_state.pop("result", None)
        st.error(f"Could not generate the DOE: {exc}")

# ---------------------------------------------------------------- results
result = st.session_state.get("result")
if result:
    runs_df = pd.DataFrame(result["runs"])
    report = result["report"]
    factor_names = list(result["factors"])

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total runs", report["total_runs"])
    m2.metric("Unique design points", report["base_design_runs"])
    m3.metric("Main-effect DOF", report["main_effect_dof"])
    m4.metric("DOF left over", report["remaining_dof_after_main_effects"])

    full_size = 1
    for q in report["factor_level_counts"].values():
        full_size *= q
    if report["method"] != "full factorial" and full_size > report["base_design_runs"]:
        st.success(
            f"{report['method'].capitalize()}: {report['base_design_runs']} runs instead of "
            f"{full_size:,} for a full factorial "
            f"({100 * report['base_design_runs'] / full_size:.2g}% of the experiments)."
        )

    tab_runs, tab_summary, tab_check, tab_download = st.tabs(
        ["Run sheet", "Design summary", "Balance check", "Download"]
    )

    with tab_runs:
        show_std = st.checkbox("Show standard run and replicate columns", value=False)
        view = runs_df.copy()
        hidden = [] if show_std else ["StandardRun", "Replicate"]
        view = view.drop(columns=[c for c in hidden if c in view.columns])
        view = view.rename(columns={"RunOrder": "Run"})
        st.dataframe(view, width="stretch", hide_index=True, height=min(700, 40 + 35 * len(view)))
        if result["response_columns"]:
            st.caption("Response columns are blank. Fill them in the downloaded file after each run.")

    with tab_summary:
        st.dataframe(report_to_frame(report), width="stretch", hide_index=True)
        st.markdown("**Orthogonal array groups**")
        if report["design_groups"]:
            groups_df = pd.DataFrame(
                [
                    {
                        "Levels": g["level_count"],
                        "Array": g["array"],
                        "Factors": ", ".join(g["factors"]),
                        "Columns used": str(g["columns_used"]),
                        "Pair frequency": g["pair_frequency"],
                    }
                    for g in report["design_groups"]
                ]
            ).astype(str)
            st.dataframe(groups_df, width="stretch", hide_index=True)
        else:
            st.write("Full factorial: every combination is run once per replicate.")
        st.markdown("**Degrees of freedom per factor**")
        st.dataframe(
            pd.DataFrame(
                {"Factor": list(report["factor_dof"]), "Main-effect DOF": list(report["factor_dof"].values())}
            ),
            width="stretch",
            hide_index=True,
        )
        st.info(
            "Leftover DOF can hold interactions or other unmodeled effects, so it isn't automatically "
            "pure error. This tool assigns main effects only. Check aliasing before estimating any "
            "interaction."
        )

    with tab_check:
        st.write("Every level of a factor should appear equally often. Counts below include replicates.")
        balance_rows = []
        for name in factor_names:
            counts = runs_df[name].value_counts().reindex(result["factors"][name])
            for level, count in counts.items():
                balance_rows.append({"Factor": name, "Level": str(level), "Runs": int(count)})
        st.dataframe(pd.DataFrame(balance_rows), width="stretch", hide_index=True)

        if len(factor_names) >= 2:
            st.write("Pairwise check: each cell should hold the same count in an orthogonal design.")
            c1, c2 = st.columns(2)
            first = c1.selectbox("Factor A", factor_names, index=0)
            second = c2.selectbox("Factor B", factor_names, index=1)
            if first == second:
                st.warning("Pick two different factors.")
            else:
                cross = pd.crosstab(runs_df[first], runs_df[second]).reindex(
                    index=result["factors"][first], columns=result["factors"][second]
                )
                st.dataframe(cross, width="stretch")
                if cross.values.min() == cross.values.max():
                    st.success(f"Balanced: every pair appears {int(cross.values.min())} time(s).")
                else:
                    st.warning("Pair counts differ. This is expected only for non-orthogonal designs.")

    with tab_download:
        out_df = runs_df.rename(columns={"RunOrder": "Run"})
        csv_bytes = out_df.to_csv(index=False).encode("utf-8-sig")
        excel_bytes = build_excel(out_df, report_to_frame(report), factor_summary(result["factors"]))
        d1, d2 = st.columns(2)
        d1.download_button(
            "Download CSV", csv_bytes, file_name="doe_table.csv", mime="text/csv", width="stretch"
        )
        d2.download_button(
            "Download Excel (run sheet + summary)",
            excel_bytes,
            file_name="doe_table.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            width="stretch",
        )
else:
    st.info("Define your factors and press **Generate DOE table** to see the run sheet.")
