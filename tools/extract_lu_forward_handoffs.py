#!/usr/bin/env python3
import argparse
import re
import sqlite3
from collections import defaultdict


FACTOR = "nvkernel_y_line_solvers_ys_factor_pipelined_batch"
FWD_FIRST = "nvkernel_y_line_solvers_ys_forward_pipelined_batch__"
FWD_CONT = "nvkernel_y_line_solvers_ys_forward_pipelined_batch_continue"


def rank_from_path(path, fallback):
    match = re.search(r"_(\d+)\.sqlite$", path)
    if match:
        return int(match.group(1))
    return fallback


def session_start_ns(con):
    try:
        row = con.execute("select utcEpochNs from TARGET_INFO_SESSION_START_TIME limit 1").fetchone()
    except sqlite3.Error:
        row = None
    return int(row[0]) if row and row[0] is not None else 0


def fetch_events(db_path, rank=None, align_session=False):
    con = sqlite3.connect(db_path)
    offset = session_start_ns(con) if align_session else 0
    rows = con.execute(
        """
        select k.start, k.end, k.globalPid, s.value
        from CUPTI_ACTIVITY_KIND_KERNEL k
        join StringIds s on k.demangledName = s.id
        where s.value like 'nvkernel_y_line_solvers_ys_factor_pipelined_batch%'
           or s.value like 'nvkernel_y_line_solvers_ys_forward_pipelined_batch%'
           or s.value like 'nvkernel_y_line_solvers_ys_forward_pipelined_batch_continue%'
        order by k.start
        """
    ).fetchall()
    con.close()
    return [
        {
            "start": int(start) + offset,
            "end": int(end) + offset,
            "pid": int(pid) if rank is None else int(rank),
            "name": name,
        }
        for start, end, pid, name in rows
    ]


def split_cases(events, gap_ns):
    cases = []
    cur = []
    last = None
    for ev in events:
        if last is not None and ev["start"] - last > gap_ns and cur:
            cases.append(cur)
            cur = []
        cur.append(ev)
        last = ev["start"]
    if cur:
        cases.append(cur)
    return cases


def rank_order_for_case(case):
    first_factor = {}
    for ev in case:
        if ev["name"].startswith(FACTOR):
            first_factor.setdefault(ev["pid"], ev["start"])
    return [pid for pid, _ in sorted(first_factor.items(), key=lambda item: item[1])]


def forward_events_by_rank(case):
    out = defaultdict(list)
    for ev in case:
        if ev["name"].startswith(FWD_FIRST) or ev["name"].startswith(FWD_CONT):
            out[ev["pid"]].append(ev)
    for pid in out:
        out[pid].sort(key=lambda ev: ev["start"])
    return out


def summarize_case(case_index, case):
    order = rank_order_for_case(case)
    forwards = forward_events_by_rank(case)
    if not order:
        return None
    nbatches = len(forwards.get(order[0], []))
    if nbatches == 0:
        return None
    handoffs = []
    for left_i in range(len(order) - 1):
        left = order[left_i]
        right = order[left_i + 1]
        left_fwd = forwards.get(left, [])
        right_fwd = forwards.get(right, [])
        for batch in range(min(len(left_fwd), len(right_fwd))):
            gap_ns = right_fwd[batch]["start"] - left_fwd[batch]["end"]
            handoffs.append(
                {
                    "from_rank": left_i,
                    "to_rank": left_i + 1,
                    "batch": batch + 1,
                    "left_end_ms": left_fwd[batch]["end"] / 1.0e6,
                    "right_start_ms": right_fwd[batch]["start"] / 1.0e6,
                    "gap_us": gap_ns / 1.0e3,
                }
            )
    return {
        "case_index": case_index,
        "nbatches": nbatches,
        "nranks": len(order),
        "rank_pids": order,
        "handoffs": handoffs,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract pipelined-LU forward handoff gaps from an Nsight Systems SQLite export."
    )
    parser.add_argument("sqlite", nargs="+")
    parser.add_argument(
        "--case-gap-ms",
        type=float,
        default=25.0,
        help="Gap between pipeline kernels that starts a new autotune candidate.",
    )
    args = parser.parse_args()

    multi_file = len(args.sqlite) > 1
    events = []
    for idx, db_path in enumerate(args.sqlite):
        rank = rank_from_path(db_path, idx) if multi_file else None
        events.extend(fetch_events(db_path, rank=rank, align_session=multi_file))
    events.sort(key=lambda ev: ev["start"])
    cases = split_cases(events, int(args.case_gap_ms * 1.0e6))
    summaries = [summarize_case(i + 1, case) for i, case in enumerate(cases)]
    summaries = [s for s in summaries if s is not None]

    print("case, batches, nranks, handoff, batch, forward_end_left_ms, forward_start_right_ms, gap_us")
    for summary in summaries:
        for h in summary["handoffs"]:
            print(
                f"{summary['case_index']}, {summary['nbatches']}, {summary['nranks']}, "
                f"{h['from_rank']}->{h['to_rank']}, {h['batch']}, "
                f"{h['left_end_ms']:.6f}, {h['right_start_ms']:.6f}, {h['gap_us']:.3f}"
            )

    by_batches = defaultdict(list)
    for summary in summaries:
        for h in summary["handoffs"]:
            by_batches[(summary["nbatches"], h["from_rank"], h["to_rank"])].append(h["gap_us"])

    print()
    print("summary_by_batches, handoff, count, min_gap_us, median_gap_us, max_gap_us")
    for (nbatches, left, right), values in sorted(by_batches.items()):
        values = sorted(values)
        median = values[len(values) // 2]
        print(
            f"{nbatches}, {left}->{right}, {len(values)}, "
            f"{values[0]:.3f}, {median:.3f}, {values[-1]:.3f}"
        )


if __name__ == "__main__":
    main()
