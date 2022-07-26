
# needs to be run from the base repo directory

# base

docker buildx build --file docker/focal/base/Dockerfile --platform linux/arm64             --tag terraferma/base:fenics-2019.1.0-focal-arm64 --push .
docker buildx build --file docker/focal/base/Dockerfile --platform linux/amd64             --tag terraferma/base:fenics-2019.1.0-focal-amd64 --push .
docker buildx build --file docker/focal/base/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/base:fenics-2019.1.0-focal       --push .

# dev-env

docker buildx build --build-arg TAG=fenics-2019.1.0-focal-arm64 --file docker/focal/dev-env/Dockerfile --platform linux/arm64             --tag terraferma/dev-env:fenics-2019.1.0-focal-arm64 --push .
docker buildx build --build-arg TAG=fenics-2019.1.0-focal-amd64 --file docker/focal/dev-env/Dockerfile --platform linux/amd64             --tag terraferma/dev-env:fenics-2019.1.0-focal-amd64 --push .
docker buildx build --build-arg TAG=fenics-2019.1.0-focal       --file docker/focal/dev-env/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/dev-env:fenics-2019.1.0-focal       --push .

# dev

docker buildx build --build-arg TAG=fenics-2019.1.0-focal-arm64 --file docker/focal/dev/Dockerfile --platform linux/arm64             --tag terraferma/dev:fenics-2019.1.0-focal-arm64 --push .
docker buildx build --build-arg TAG=fenics-2019.1.0-focal-amd64 --file docker/focal/dev/Dockerfile --platform linux/amd64             --tag terraferma/dev:fenics-2019.1.0-focal-amd64 --push .
docker buildx build --build-arg TAG=fenics-2019.1.0-focal       --file docker/focal/dev/Dockerfile --platform linux/amd64,linux/arm64 --tag terraferma/dev:fenics-2019.1.0-focal       --push .

