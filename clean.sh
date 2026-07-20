echo ${BASH_SOURCE}
BASEDIR=$(dirname "$0")
#echo "$BASEDIR"
cd "$BASEDIR/src" || exit

find . -type f -name '*.o' -delete
find ./ -type f -name '*.so' -delete
find ./ -type f -name '*.a' -delete
