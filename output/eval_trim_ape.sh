#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# 用法说明
# ==========================================================
usage() {
cat << EOF
Usage:
  $0 EST_TRAJ GT_TRAJ [CUT_FRONT] [CUT_BACK]

Arguments:
  EST_TRAJ    Estimated trajectory (TUM format or superset)
  GT_TRAJ     Ground-truth trajectory (TUM format)
  CUT_FRONT   Seconds to trim at the beginning (default: 3)
  CUT_BACK    Seconds to trim at the end (default: 2)

Example:
  $0 pose.txt MH_01.txt 3 2
EOF
exit 1
}

# ==========================================================
# 参数解析
# ==========================================================
[[ $# -lt 2 ]] && usage

EST_TRAJ="$1"
GT_TRAJ="$2"
CUT_FRONT="${3:-3}"
CUT_BACK="${4:-2}"

OUT_EST="est_trim.txt"
OUT_GT="gt_trim.txt"

# ==========================================================
# 工具函数
# ==========================================================
first_ts() { awk 'NF>=1 && $1 ~ /^[0-9.eE+-]+$/ {print $1; exit}' "$1"; }
last_ts()  { awk 'NF>=1 && $1 ~ /^[0-9.eE+-]+$/ {ts=$1} END{print ts}' "$1"; }

sanitize_tum8() {
  local in="$1"
  local out="$2"
  awk '
    BEGIN{OFS=" "}
    NF>=8 {
      print $1,$2,$3,$4,$5,$6,$7,$8
    }
  ' "$in" | sed -E 's/[[:space:]]+$//' > "$out"
}

# ==========================================================
# Step 1: 清洗输入为严格 TUM 8 列
# ==========================================================
TMP_EST=".tmp_est_tum8.txt"
TMP_GT=".tmp_gt_tum8.txt"

sanitize_tum8 "$EST_TRAJ" "$TMP_EST"
sanitize_tum8 "$GT_TRAJ"  "$TMP_GT"

# ==========================================================
# Step 2: 计算时间窗口（交集 + 前后裁剪）
# ==========================================================
EST_START=$(first_ts "$TMP_EST")
EST_END=$(last_ts  "$TMP_EST")
GT_START=$(first_ts "$TMP_GT")
GT_END=$(last_ts  "$TMP_GT")

START_TIME=$(awk -v a="$EST_START" -v b="$GT_START" 'BEGIN{print (a>b)?a:b}')
END_TIME=$(awk -v a="$EST_END"   -v b="$GT_END"   'BEGIN{print (a<b)?a:b}')

START_TIME=$(awk -v s="$START_TIME" -v c="$CUT_FRONT" 'BEGIN{printf "%.9f\n", s+c}')
END_TIME=$(awk -v e="$END_TIME"   -v c="$CUT_BACK"  'BEGIN{printf "%.9f\n", e-c}')

echo "=================================================="
echo "[INFO] EST raw range : $EST_START  ~  $EST_END"
echo "[INFO]  GT raw range : $GT_START   ~  $GT_END"
echo "[INFO] Trim settings : front=${CUT_FRONT}s, back=${CUT_BACK}s"
echo "[INFO] Eval window   : $START_TIME  ~  $END_TIME"
echo "=================================================="

awk -v s="$START_TIME" -v e="$END_TIME" 'BEGIN{ if (!(e>s)) exit 1 }' || {
  echo "[FATAL] Invalid evaluation window."
  echo "        Check timestamp unit / time base of EST vs GT."
  exit 2
}

# ==========================================================
# Step 3: 裁剪轨迹
# ==========================================================
awk -v s="$START_TIME" -v e="$END_TIME" \
    '($1>=s && $1<=e){print}' "$TMP_EST" > "$OUT_EST"

awk -v s="$START_TIME" -v e="$END_TIME" \
    '($1>=s && $1<=e){print}' "$TMP_GT" > "$OUT_GT"

N_EST=$(wc -l < "$OUT_EST")
N_GT=$(wc -l < "$OUT_GT")

echo "[INFO] Trimmed rows: est=$N_EST, gt=$N_GT"

if [[ $N_EST -lt 10 || $N_GT -lt 10 ]]; then
  echo "[WARN] Very few poses after trimming. APE may be unreliable."
fi

# ==========================================================
# Step 4: evo APE
# ==========================================================
evo_ape tum "$OUT_GT" "$OUT_EST" \
    -a \
    --t_max_diff 0.01 \
    -v --plot

