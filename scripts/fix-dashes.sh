#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Convert hyphens between capitalized names (e.g, Gram-Schmidt) into en dashes
# (i.e., Gram--Schmidt) in .rst, .tex and .bib files.
#
# Usage:
#
# ./fix-dashes.sh terms   [paths] list each term that would change, with counts
# ./fix-dashes.sh preview [paths] show changes with file and line number
# ./fix-dashes.sh apply   [paths] edit files in place
#
# Paths defaults to the current directory.
# ------------------------------------------------------------------------------

# Exit on errors, unset variables, and failures anywhere in a pipeline.
set -euo pipefail

# --- Requirements ------------------------------------------------------------

# rg finds the files and previews changes.
if ! command -v rg >/dev/null 2>&1; then
  echo "Error: ripgrep (rg) is not installed or not in PATH." >&2
  echo "Install it from https://github.com/BurntSushi/ripgrep#installation" >&2
  exit 2
fi

if ! rg --pcre2-version >/dev/null 2>&1; then
  echo "Error: rg was built without PCRE2 support." >&2
  echo "Install a build with PCRE2." >&2
  exit 2
fi

# perl does the in-place edits.
if ! command -v perl >/dev/null 2>&1; then
  echo "Error: perl is not installed or not on PATH." >&2
  exit 2
fi

# --- Skip lists: edit these -------------------------------------------------

# Prefix words: a hyphen directly after one of these is never converted
SKIP_WORDS=(
  Band
  Bi
  Block
  Computer
  Differential
  Distributed
  Edge
  Entropy
  Fixed
  Higher
  Initial
  Large
  Least
  Low
  Matrix
  Multi
  Next
  Non
  Post
  Pre
  Quasi
  Reduced
  Right
  Sensitivity
  Stability
  Step
  Structure
  Three
  Time
  Transpose
  Type
  User
  Variable
)

# Exact word pairs left alone: The second word must match exactly. End it with *
# to match any word starting with it instead. LaTeX accents can be written
# as-is: the second word is cut at its first non-letter
SKIP_PAIRS=(
  "Advection-Dominated"
  "Almgren-Bell"
  "Berger-Vergiat"
  "Chi-Wang"
  "Frank-Kam*"
  "Lebrun-Grandi*"
  "Peer-Timo"
  "Pollaczek-Geiringer"
  "Quintana-Ort*"
  "Sacks-Davis"
  "Springer-Verlag"
  "Swig-Fortran"
  "Trilinos-Overview"
)

# File types to process.
GLOBS=( '*.rst' '*.tex' '*.bib' )

# --- Build the pattern -------------------------------------------------------

# Anything inside a web address is skipped entirely. The URL stops at
# whitespace, braces, angle brackets or quotes, so \href{url}{Text} and `text
# <url>`_ still convert the link text.
#
# Then: a capitalized word (not the "License" in SPDX-License-Identifier),
# which \K drops from the match so that chains like A-B-C work.
PAT='https?://[^\s{}<>"]+(*SKIP)(*FAIL)|(?<!SPDX-)\b[A-Z][a-z]+\K'

# Skip the hyphen if the word before it is exactly a SKIP_WORDS entry.
for w in "${SKIP_WORDS[@]}"; do
  PAT+="(?<!\\b${w})"
done

# Skip the hyphen if the words around it form a SKIP_PAIRS entry. A trailing *
# drops the word-end check so the second word acts as a prefix.
for p in "${SKIP_PAIRS[@]}"; do
  first=${p%%-*}
  second=${p#*-}
  if [[ $second == *'*' ]]; then
    second=${second%'*'}
    end=''
  else
    end='\b'
  fi
  second=${second%%[!A-Za-z]*}
  PAT+="(?!(?<=\\b${first})-${second}${end})"
done

# The hyphen itself, followed by another capitalized word.
PAT+='-(?=[A-Z][a-z]+\b)'

# Exported so perl can read it as $ENV{PAT} without shell quoting issues.
export PAT

# --- Run ---------------------------------------------------------------------

# First argument is the mode; the rest are paths (default: current directory).
mode=${1:-}
shift || true
paths=( "$@" )
[ ${#paths[@]} -eq 0 ] && paths=( . )

# Unknown or missing mode: print the usage block from the header and stop.
case "$mode" in
  terms|preview|apply) ;;
  *) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac

# Turn GLOBS into rg arguments: -g '*.rst' -g '*.tex' ...
glob_args=()
for g in "${GLOBS[@]}"; do glob_args+=( -g "$g" ); done

# Stop early if nothing matches, so xargs never runs perl with no files.
if ! rg -P -q "${glob_args[@]}" "$PAT" "${paths[@]}"; then
  echo "No matches." >&2
  exit 0
fi

case "$mode" in
  terms)
    # For each line that would change, print the converted terms with the
    # occurrence count and ordered by frequency.
    #
    # The $ in the perl code is Perl syntax and must not expand in the shell.
    # shellcheck disable=SC2016
    rg -P -l0 "${glob_args[@]}" "$PAT" "${paths[@]}" \
      | xargs -0 perl -ne '
          my $o = $_; s/$ENV{PAT}/--/g; next if $_ eq $o;
          print "$1\n" while /(\b[A-Z][a-z]+(?:-{1,2}[A-Z][a-z]+)*--[A-Z][a-z]+(?:-{1,2}[A-Z][a-z]+)*\b)/g;
        ' \
      | sort | uniq -c | sort -rn
    ;;
  preview)
    # Show changes with file and line number; files are not modified.
    rg -P "${glob_args[@]}" "$PAT" --replace='--' "${paths[@]}"
    ;;
  apply)
    # rg lists matching files (NUL-separated); perl edits them in place.
    #
    # Each cell in an .rst table is handled separately to maintain | alignment:
    # count the replacements, then remove that many padding spaces before the
    # closing |, always keeping at least one. If there is not enough padding,
    # warn with file and line number. "close ARGV if eof" resets the line
    # counter $. for each file so warnings report correct line numbers.
    #
    # The $ in the perl code is Perl syntax and must not expand in the shell.
    # shellcheck disable=SC2016
    rg -P -l0 "${glob_args[@]}" "$PAT" "${paths[@]}" \
      | xargs -0 perl -pi -e '
          if ($ARGV =~ /\.rst$/ && /\|\s*$/) {
            s{([^|]*)\|}{
              my $c = $1;
              my $n = ($c =~ s/$ENV{PAT}/--/g) || 0;
              if ($n) {
                if ($c =~ / {$n} $/) { $c =~ s/ {$n}( )$/$1/ }
                else { warn "$ARGV:$.: not enough padding, table may be misaligned\n" }
              }
              "$c|"
            }ge;
          } else {
            s/$ENV{PAT}/--/g;
          }
          close ARGV if eof;
        '
    ;;
esac
