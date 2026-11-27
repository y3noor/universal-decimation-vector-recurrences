#!/bin/sh
set -eu
cd "$(dirname "$0")"
case "$#:$*" in
    "1:records")
        exec python3 verify_searches.py --records-only
        ;;
    "3:case "*)
        if [ "$3" = 2 ]; then
            program=search_quadratic
        else
            program=search_general
        fi
        c++ -O3 -std=c++17 "$program.cpp" -o "$program"
        exec python3 verify_searches.py --case "$2" "$3"
        ;;
    "1:quick"|"1:full") mode=$1 ;;
    *) printf 'Usage: sh replay.sh quick|full|records|case N M\n' >&2; exit 2 ;;
esac
for program in search_quadratic search_general; do
    c++ -O3 -std=c++17 "$program.cpp" -o "$program"
done
python3 verify_feedback.py
python3 verify_independent.py
python3 verify_sizes.py
python3 verify_cubic_example.py
if [ "$mode" = quick ]; then
    python3 verify_searches.py --quick
else
    python3 verify_searches.py
fi
