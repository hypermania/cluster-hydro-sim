#!/usr/bin/env python3
"""Run all eight Heggie gas models with bounded single-threaded concurrency."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--executable", type=Path, default=Path("main-strict"))
    parser.add_argument("--zones", type=int, default=400)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.jobs <= 6 or args.zones < 3:
        parser.error("require 1..6 jobs and at least 3 zones")
    args.output.mkdir(parents=True, exist_ok=True)
    cases = [(f"single_{f:g}", 0, f) for f in (.001, .002, .005, .01, .02)]
    cases += [("segregation", 1, .005), ("bs", 2, .005), ("bsbb", 3, .005)]
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")

    def run(case):
        name, model, fstar = case
        directory = args.output/name
        # Do not overwrite existing histories, including failed/partial runs.
        if directory.exists():
            raise FileExistsError(directory)
        with (args.output/f"{name}.log").open("x") as log:
            code = subprocess.call([str(args.executable.resolve()), "moving-heggie",
                                    str(directory), str(model), str(fstar), str(args.zones)],
                                   stdout=log, stderr=subprocess.STDOUT, env=env)
        print(f"{name}: exit={code}", flush=True)
        return code

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        codes = list(pool.map(run, cases))
    raise SystemExit(int(any(codes)))


if __name__ == "__main__":
    main()
