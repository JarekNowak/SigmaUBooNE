#!/bin/bash
# In-container part of fetch_sample.sh (SL7 container: /usr/bin/xrdcp and xrdfs v5, which support
# bearer tokens; the larsoft xrootd 4.8 does not). Not run directly. Args: DEST SUBDEF RAWDIMS MAX PAR
#   SUBDEF/RAWDIMS non-empty -> (re)build DEST/pnfs.list from SAM: the sub-definition if it exists,
#   else the raw dimensions; MAX -1 = all online, 0 = list only
DEST="$1"; SUBDEF="$2"; RAW="$3"; MAX="$4"; PAR="$5"
DOOR=root://fndca1.fnal.gov:1094
[ -z "$BEARER_TOKEN_FILE" ] && unset BEARER_TOKEN_FILE
WORK=$(mktemp -d "$DEST/.work_XXXX"); trap 'rm -rf "$WORK"' EXIT

# 1. SAM -> dCache paths (prefer the dCache copy that prestage-dataset makes over the tape path)
if [ -n "$SUBDEF$RAW" ]; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh >/dev/null 2>&1
  setup sam_web_client >/dev/null 2>&1 || { echo "ERROR: cannot set up sam_web_client" >&2; exit 1; }
  if [ -n "$SUBDEF" ] && samweb -e uboone describe-definition "$SUBDEF" >/dev/null 2>&1; then DIMS="defname: $SUBDEF"
  else DIMS="$RAW"; echo "  sub-definition $SUBDEF not found (not prestaged?): using the sample dimensions"; fi
  samweb -e uboone list-files "$DIMS" > "$WORK/names" || exit 1
  echo "  SAM: $(wc -l < "$WORK/names") files in '$DIMS'; locating (8 at a time) ..."
  xargs -P 8 -I{} sh -c '
    locs=$(samweb -e uboone locate-file "$1")
    d=$(echo "$locs" | grep -o "^dcache:/pnfs/[^(]*" | head -1 | sed "s/^dcache://")
    [ -n "$d" ] || d=$(echo "$locs" | grep -o "/pnfs/[^(]*" | head -1)
    [ -n "$d" ] && echo "$d/$1"' _ {} < "$WORK/names" | sort > "$DEST/pnfs.list"
fi
[ -s "$DEST/pnfs.list" ] || { echo "ERROR: empty $DEST/pnfs.list" >&2; exit 1; }

# 2. not yet local
while read -r p; do [ -s "$DEST/$(basename "$p")" ] || echo "$p"; done < "$DEST/pnfs.list" > "$WORK/todo"
nlocal=$(( $(wc -l < "$DEST/pnfs.list") - $(wc -l < "$WORK/todo") ))

# 3. which are on disk at FNAL (a tape path whose dCache copy is online is switched to it);
#    ERR = stat failed (expired token, no permission, door down): not the same as tape-only
xargs -P 16 -I{} sh -c '
  p="$1"; st=ERR
  for q in "$p" "$(echo "$p" | sed "s|^/pnfs/uboone/data/|/pnfs/uboone/production/data/|")"; do
    f=$(/usr/bin/xrdfs '"$DOOR"' stat "/pnfs/fnal.gov/usr/${q#/pnfs/}" 2>/dev/null | grep "^Flags:")
    [ -n "$f" ] && st=OFF
    if [ -n "$f" ] && [ "${f#*Offline}" = "$f" ]; then echo "ON $q"; exit 0; fi
    [ "$q" = "$p" ] || break
  done
  echo "$st $p"' _ {} < "$WORK/todo" > "$WORK/state"
nerr=$(grep -c '^ERR ' "$WORK/state")
[ "$nerr" -gt 0 ] && echo "  WARNING: stat failed for $nerr files (expired token? run get_token.sh)" >&2
grep '^ON ' "$WORK/state" | cut -c4- > "$WORK/copy"
noff=$(grep -c '^OFF ' "$WORK/state")
[ "$MAX" -ge 0 ] && head -n "$MAX" "$WORK/copy" > "$WORK/c" && mv "$WORK/c" "$WORK/copy"
echo "  local $nlocal   online at FNAL $(grep -c '^ON ' "$WORK/state")   tape-only $noff   stat error $nerr   copying $(wc -l < "$WORK/copy") ($PAR parallel)"
[ "$MAX" = 0 ] && exit 0

# 4. copy to .part, rename on success; one retry; errors to DEST/.xrdcp_errors.log
xargs -P "$PAR" -I{} sh -c '
  p="$1"; dst="'"$DEST"'/$(basename "$p")"
  [ -s "$dst" ] && exit 0
  for try in 1 2; do
    err=$(/usr/bin/xrdcp --nopbar --retry 2 -f "'"$DOOR"'/pnfs/fnal.gov/usr/${p#/pnfs/}" "$dst.part" 2>&1 >/dev/null) \
      && mv "$dst.part" "$dst" && { echo "ok   $(basename "$p")"; exit 0; }
    rm -f "$dst.part"
    echo "$(date "+%F %T") try $try $(basename "$p"): $(echo "$err" | tail -2 | tr "\n" " ")" >> "'"$DEST"'/.xrdcp_errors.log"
    [ $try -eq 1 ] && sleep 30
  done
  echo "FAIL $p" >&2' _ {} < "$WORK/copy" > "$WORK/copied"
echo "  copied $(grep -c '^ok ' "$WORK/copied")   failed $(( $(wc -l < "$WORK/copy") - $(grep -c '^ok ' "$WORK/copied") ))   tape-only (prestage first) $noff"
