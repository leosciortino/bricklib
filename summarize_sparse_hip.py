#!/usr/bin/env python3
"""Summarize sparse BrickLib rocprofv3 counters without plotting."""

import argparse
import csv
from pathlib import Path


COUNTERS = {
    "FETCH_SIZE",
    "WRITE_SIZE",
    "SQ_INSTS_VALU_ADD_F64",
    "SQ_INSTS_VALU_MUL_F64",
    "SQ_INSTS_VALU_FMA_F64",
    "SQ_INSTS_VALU_TRANS_F64",
}


def summarize(path):
    dispatches = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if "_brick_trans" not in row["Kernel_Name"]:
                continue
            dispatch = dispatches.setdefault(
                int(row["Dispatch_Id"]),
                {
                    "start": float(row["Start_Timestamp"]),
                    "end": float(row["End_Timestamp"]),
                },
            )
            dispatch[row["Counter_Name"]] = float(row["Counter_Value"])

    runs = [dispatches[key] for key in sorted(dispatches)]
    if len(runs) != 101:
        raise ValueError(f"{path}: expected 1 warmup + 100 runs, found {len(runs)}")
    runs = runs[1:]

    for run in runs:
        missing = COUNTERS - run.keys()
        if missing:
            raise ValueError(f"{path}: missing counters: {', '.join(sorted(missing))}")

    total_bytes = 1024 * sum(run["FETCH_SIZE"] + run["WRITE_SIZE"] for run in runs)
    total_time = sum(run["end"] - run["start"] for run in runs) * 1e-9
    total_flop = 64 * sum(
        run["SQ_INSTS_VALU_ADD_F64"]
        + run["SQ_INSTS_VALU_MUL_F64"]
        + 2 * run["SQ_INSTS_VALU_FMA_F64"]
        + run["SQ_INSTS_VALU_TRANS_F64"]
        for run in runs
    )
    return total_bytes, total_time, total_flop


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results",
        nargs="?",
        type=Path,
        default=Path("/p/vast1/sciortino2/bricklib_line_cube"),
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("bricklib_sparse.csv"))
    args = parser.parse_args()

    records = []
    for workload in sorted(path for path in args.results.iterdir() if path.is_dir()):
        profiles = list(workload.rglob("*counter_collection.csv"))
        if not profiles:
            continue
        if len(profiles) != 1:
            raise ValueError(f"{workload}: expected one counter CSV, found {len(profiles)}")
        total_bytes, total_time, total_flop = summarize(profiles[0])
        records.append(
            ["bricklib", workload.name, total_bytes, total_time, total_flop]
        )

    if not records:
        raise ValueError(f"no rocprofv3 counter CSV files found below {args.results}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["library", "workload_name", "total_bytes", "total_time", "total_flop"])
        writer.writerows(records)
    print(f"Wrote {len(records)} workloads to {args.output}")


if __name__ == "__main__":
    main()
