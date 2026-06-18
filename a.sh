FILES=( $(ls runs/N*/run_*/output.txt | sort -V) )
NS=$(printf "%s\n" "${FILES[@]}" \
       | sed -E 's|.*/N([0-9]+)/run_.*|\1|' \
       | paste -sd, -)
python3 python/analyze.py "${FILES[@]}" \
    --Ns "$NS" \
    --r_outer 40 --r_inner 1 --radius 1 --dS 0.2 \
    --out outputs/tp3_plots/sistema1