# Script to build (and push) docker images for ThermoCodegen.

# Note that this assumes that an appropriate multiplatform docker builder is available.
# i.e. it may be necessary to run:
# > docker buildx create --use
# first.

usage() {
    echo "Usage:" 1>&2
    echo "bash $0 [-h] [-t tag<string>] [-b branch<string>] [-r repo<string>] [-d] [-a] dir<string>" 1>&2
    echo "  dir: required name of the subdirectory containing the Dockerfile" 1>&2
    echo "  -h: print this help message and exit" 1>&2
    echo "  -t: specify a tag name (defaults to fenics-2019.1.0-focal)" 1>&2
    echo "  -b: specify a branch name (defaults to master-2019.1.0)" 1>&2
    echo "  -r: specify a repo URL (defaults to https://github.com/TerraFERMA/TerraFERMA.git)" 1>&2
    echo "  -d: enable debugging (if tag name is default, suffixes name with -debug)" 1>&2
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
BRANCH='master-2019.1.0'
REPO='https://github.com/TerraFERMA/TerraFERMA.git'
DEBUG=0
PLATFORMS=0

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
           PLATFORMS=${OPTARG}
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

# if no tag is specified default to enki-portal
if [ -z "$TAG" ]; then
    TAG='fenics-2019.1.0-focal'
    if [ "$BRANCH" != 'master-2019.1.0' ]; then
        TAG="${TAG}-$BRANCH"
    fi
    if [ "$BUILD" == 'Debug' ]; then
        TAG="${TAG}-debug"
    fi
fi

echo $DIR
echo $SDIR
echo $DEBUG
echo $TAG
echo $REPO

cd $repo_path
docker buildx build --file $DIR/Dockerfile \
                    --build-arg TAG=$TAG \
                    --build-arg BRANCH=$BRANCH \
                    --build-arg DEBUG=$DEBUG \
                    --build-arg REPO=$REPO \
                    --tag terraferma/$SDIR:$TAG --push .


## needs to be run from the base repo directory
#
## base
#
#docker buildx build --file docker/focal/base/Dockerfile --platform linux/arm64             --tag terraferma/base:fenics-2019.1.0-focal-arm64 --push .
#docker buildx build --file docker/focal/base/Dockerfile --platform linux/amd64             --tag terraferma/base:fenics-2019.1.0-focal-amd64 --push .
#docker buildx build --file docker/focal/base/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/base:fenics-2019.1.0-focal       --push .
#
## dev-env
#
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal-arm64 --file docker/focal/dev-env/Dockerfile --platform linux/arm64             --tag terraferma/dev-env:fenics-2019.1.0-focal-arm64 --push .
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal-amd64 --file docker/focal/dev-env/Dockerfile --platform linux/amd64             --tag terraferma/dev-env:fenics-2019.1.0-focal-amd64 --push .
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal       --file docker/focal/dev-env/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/dev-env:fenics-2019.1.0-focal       --push .
#
## dev
#
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal-arm64 --file docker/focal/dev/Dockerfile --platform linux/arm64             --tag terraferma/dev:fenics-2019.1.0-focal-arm64 --push .
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal-amd64 --file docker/focal/dev/Dockerfile --platform linux/amd64             --tag terraferma/dev:fenics-2019.1.0-focal-amd64 --push .
#docker buildx build --build-arg TAG=fenics-2019.1.0-focal       --file docker/focal/dev/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/dev:fenics-2019.1.0-focal       --push .
#
