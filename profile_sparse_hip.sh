#!/usr/bin/env bash

# Profile an already-built sparse HIP executable on one GPU.
set -euo pipefail

EXE="${EXE:-./build-hip/single/sparse-hip}"
OUT_DIR="${OUT_DIR:-/p/vast1/sciortino2/bricklib_line_cube}"
WIDTH="${WIDTH:-1}"

PASSES=(
    "fetch FETCH_SIZE"
    "write WRITE_SIZE"
    "flop SQ_INSTS_VALU_ADD_F64 SQ_INSTS_VALU_MUL_F64 SQ_INSTS_VALU_FMA_F64 SQ_INSTS_VALU_TRANS_F64"
)

command -v rocprofv3 >/dev/null || {
    echo "rocprofv3 was not found" >&2
    exit 1
}
[[ -x "$EXE" ]] || {
    echo "HIP executable was not found: $EXE" >&2
    exit 1
}

mkdir -p "$OUT_DIR"

for config in "cube 1" "cube 2" "star 1" "star 2" "star 3" "star 4"; do
    read -r stencil radius <<< "$config"
    for shape in line circle; do
        run_dir="${OUT_DIR}/${stencil}_r${radius}_${shape}"
        echo "Profiling $stencil radius $radius, $shape -> $run_dir"
        for pass in "${PASSES[@]}"; do
            read -ra args <<< "$pass"
            pass_dir="${run_dir}/${args[0]}"
            mkdir -p "$pass_dir"
            rocprofv3 \
                --pmc "${args[@]:1}" \
                --output-format csv \
                --output-directory "$pass_dir" \
                --output-file profile \
                -- "$EXE" "$stencil" "$radius" "$WIDTH" "$shape"
        done
    done
done
