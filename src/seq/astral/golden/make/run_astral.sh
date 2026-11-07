#!/bin/bash
# Runs the real ASTRAL 5.7.8 (a black box, only to compare against: its code is never read or copied) on generated data
# sets and keeps what it prints, for assemble.py.
# usage: ASTRAL_JAR=/path/to/astral.5.7.8.jar run_astral.sh DATADIR REFDIR NAME:EXACT(0|1)...
# DATADIR/NAME has genes.tre, true.tre and random.tre (made by gen.py and randsp.py); REFDIR gets one folder per NAME.
# The jar is a download from the ASTRAL project (Apache-2.0): keep it in a folder of its own and run this script from
# another one, as for any downloaded program.
A="${ASTRAL_JAR:?set ASTRAL_JAR to the path of astral.5.7.8.jar}"
data="$1"; ref="$2"; shift 2
for spec in "$@"; do
  d="${spec%%:*}"; exact="${spec##*:}"
  mkdir -p "$ref/$d"
  g="$data/$d/genes.tre"
  java -jar "$A" -i "$g" -o "$ref/$d/default.tre" 2> "$ref/$d/default.log"
  java -jar "$A" -i "$g" -p 0 -o "$ref/$d/p0.tre" 2> "$ref/$d/p0.log"
  if [ "$exact" = "1" ]; then
    java -jar "$A" -i "$g" -x -o "$ref/$d/x.tre" 2> "$ref/$d/x.log"
  fi
  java -jar "$A" -q "$data/$d/true.tre" -i "$g" -t 2 -o "$ref/$d/q_true.tre" 2> "$ref/$d/q_true.log"
  k=0
  while read -r line; do
    k=$((k+1))
    echo "$line" > "$ref/$d/sp$k.tre"
    java -jar "$A" -q "$ref/$d/sp$k.tre" -i "$g" -t 2 -o "$ref/$d/q_sp$k.tre" 2> "$ref/$d/q_sp$k.log"
  done < "$data/$d/random.tre"
  def=$(grep 'Final quartet score is' "$ref/$d/default.log" | head -1 | sed 's/.*is: //')
  p0=$(grep 'Final quartet score is' "$ref/$d/p0.log" | head -1 | sed 's/.*is: //')
  ex=""
  if [ "$exact" = "1" ]; then ex=$(grep 'Final quartet score is' "$ref/$d/x.log" | head -1 | sed 's/.*is: //'); fi
  tr=$(grep 'Final quartet score is' "$ref/$d/q_true.log" | head -1 | sed 's/.*is: //')
  norm=$(grep 'Final normalized' "$ref/$d/default.log" | head -1 | sed 's/.*: //')
  echo "$d default=$def p0=$p0 exact=$ex true=$tr norm=$norm"
done
