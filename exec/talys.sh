#!/bin/sh

##################################################
#
#  cleans up after a finished TALYS calculation
#  has to be executed from within the result directory
#
##################################################

# stop the script if not in a calculation diretory
if [ ! -f calc.status ]; then
  exit 1
fi
# or if an archive already exists
if [ -f archive.tgz ]; then
  exit 1
fi

# extract head and tail from output
tail -n3 output > output_tail
awk '/### RESULTS FOR E=/{print NR-1; exit 0;}' output | xargs -I NR head -nNR output > output_head

# trim the output file
cat output_head output_tail > output
rm output_head output_tail

# keep a list of files
find -type f -printf '%P\n' | sort > filelist

# some files also remain outside archive
tmpdir=`mktemp -d`
cp input output energies filelist calc.status "$tmpdir"

# make the archive
tar -zcf archive.tgz --files-from "$tmpdir/filelist" --remove-files

# copy back
cp $tmpdir/* .

# clean up
rm -r "$tmpdir"

# check status of calculation
if grep -q "TALYS team congratulates" output; then
  echo "success" > calc.status
else
  echo "failure" > calc.status
fi



