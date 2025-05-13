# Script to build (and push) docker images for TerraFERMA.

# Note that this assumes that an appropriate multiplatform docker builder is available.
# i.e. it may be necessary to run:
# > docker buildx create --use
# first.

usage() {
    echo "Usage:" 1>&2
    echo "bash $0 [-h] [-t tag<string>] [-b branch<string>] [-r repo<string>] [-p platform<string>,platform<string>] [-d] dir<string>" 1>&2
    echo "  dir: required name of the subdirectory containing the Dockerfile" 1>&2
    echo "  -h: print this help message and exit" 1>&2
    echo "  -t: specify a tag name (defaults to <os>-<branch if not main>-<debug>-<platform>)" 1>&2
    echo "  -b: specify a branch name (defaults to main)" 1>&2
    echo "  -r: specify a repo URL (defaults to https://github.com/TerraFERMA/TerraFERMA.git)" 1>&2
    echo "  -d: enable debugging (if tag name is default, suffixes tag with -debug)" 1>&2
    echo "  -p: comma separated list of platforms (defaults to current platform)" 1>&2
}

error() {
    usage
    exit 1
}

# realpath not available by default on macs so define it here
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

full_path=$(realpath $0)
script_path=$(dirname $full_path)
repo_path=$(dirname $script_path)

# parse the arguments
TAG=''
BRANCH='main'
REPO='https://github.com/TerraFERMA/TerraFERMA.git'
DEBUG=0
PLATFORMS=''
BUILDER=''

while getopts ":t:b:r:p:dh" opt; do
    case $opt in
        h )
           usage
           exit 0
           ;;
        t )
           TAG=${OPTARG}
           ;;
        b )
           BRANCH=${OPTARG}
           ;;
        r )
           REPO=${OPTARG}
           ;;
        d )
           DEBUG=1
           ;;
        p )
           PLATFORMS="--platform ${OPTARG}"
           ;;
        : )
           echo "ERROR: -${OPTARG} requires an argument." 1>&2
           error
           ;;
        * )
           echo "ERROR: unknown option -${OPTARG}." 1>&2
           error
           ;;
    esac
done

shift $((OPTIND - 1))

if [ $# == 0 ]; then
    echo "ERROR: missing required subdirectory." 1>&2
    error
fi
DIR=$(realpath $1)
SDIR=$(basename $DIR)
ADIR=$(basename $(dirname $DIR) )

PTAG='userplatform'
if [ -z "$PLATFORMS" ]; then
    PROC=`uname -m`
    if [ "$PROC" == "x86_64" ]; then
        PTAG="amd64"
    elif [ "$PROC" == "arm64" ]; then
        PTAG="arm64"
    fi
else
    BUILDER='buildx'
fi

# if no tag is specified default to arch dir
if [ -z "$TAG" ]; then
    TAG="$ADIR"
    if [ "$BRANCH" != 'main' ]; then
        TAG="${TAG}-$BRANCH"
    fi
    if [ $DEBUG -eq 1 ]; then
        TAG="${TAG}-debug"
    fi
    if [ "$PTAG" ]; then
        TAG="${TAG}-${PTAG}"
    fi
fi

echo "Building:"
echo "  $DIR/Dockerfile"
echo "  with tag ghcr.io/terraferma/$SDIR:$TAG"
echo "  from repo $REPO"
if [ "$PLATFORMS" ]; then
  echo "  with $PLATFORMS"
fi
if [ $DEBUG -eq 1 ]; then
  echo "  with debugging enabled"
fi

cd $repo_path
docker $BUILDER build --file $DIR/Dockerfile \
                    --build-arg TAG=$TAG \
                    --build-arg BRANCH=$BRANCH \
                    --build-arg DEBUG=$DEBUG \
                    --build-arg REPO=$REPO \
                    --tag ghcr.io/terraferma/$SDIR:$TAG $PLATFORMS --push .

