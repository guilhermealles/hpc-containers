#!/bin/bash
DIR=$(dirname $0)

function usage()
{
    echo "$0 <sha1_before> [sha1_now]";
    echo "  where <sha1_before> is the SHA1 before";
    echo "  where [sha1_now] (optional) is the current SHA1 (or master)";
}

##############################
# Parameter handling section #
##############################
SHA1B=${1:-}
if [ -z "$SHA1B" ]; then
    echo "Error: <sha1_before> is empty"
    usage;
    exit;
fi

SHA1N=${2:-}
if [ -z "$SHA1N" ]; then
    SHA1N="master"
fi

function export()
{
    INPUT="2018_WSCAD.org"
    emacs -batch -l ~/.emacs.d/init.el \
	  --eval "(setq enable-local-eval t)" \
	  --eval "(setq enable-local-variables t)" \
	  --eval "(setq org-export-babel-evaluate t)" \
	  $INPUT \
	  --funcall org-latex-export-to-latex
    OUTPUT=$(basename $INPUT .org)_$1.tex
    mv $(basename $INPUT .org).tex $OUTPUT
    echo $OUTPUT
}

git checkout $SHA1B
BEFORE=$(export "before")
git checkout $SHA1N
NOW=$(export "now")

latexdiff $BEFORE $NOW > diff.tex
rubber --pdf diff.tex
