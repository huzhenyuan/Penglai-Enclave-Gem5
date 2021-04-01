#!/bin/bash

function print_usage() {
	RED='\033[0;31m'
	BLUE='\033[0;34m'
	BOLD='\033[1m'
	NONE='\033[0m'

	echo -e "\n${RED}Usage${NONE}:
	.${BOLD}/docker_cmd.sh${NONE} [OPTION]"

	echo -e "\n${RED}OPTIONS${NONE}:
	${BLUE}build${NONE}: build penglai-demo image
	${BLUE}run-qemu${NONE}: run penglai-demo image in (modified) qemu 
	"
}

if [[ $1 == *"help"* ]]; then
	print_usage
	exit 0
fi

# build gem5 
if [[ $1 == "build" ]]; then
	echo "Build: building penglai demo image"
	docker run -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.5 bash scripts/build.sh
	exit 0
fi

# run test 
if [[ $1 == "run" ]]; then
	echo "Run: run penglai demo image in sPMP-supported Qemu"
	docker run -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.5 bash scripts/run.sh
	exit 0
fi

# show the original result 
if [[ $1 == *"show"* ]]; then
	echo "Run: run docker"
	#sudo docker run --privileged --cap-add=ALL  -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.2
	docker run -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.5 bash scripts/show.sh
	exit 0
fi

# dmoe: hello-world
if [[ $1 == *"hello"* ]]; then
	echo "Run: run docker"
	#sudo docker run --privileged --cap-add=ALL  -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.2
	docker run -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.5 bash scripts/hello.sh
	exit 0
fi

# run docker 
if [[ $1 == *"docker"* ]]; then
	echo "Run: run docker"
	#sudo docker run --privileged --cap-add=ALL  -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.2
	docker run -v $(pwd):/home/penglai/penglai-enclave -w /home/penglai/penglai-enclave --rm -it ddnirvana/penglai-enclave:v0.5 bash scripts/source.sh
	exit 0
fi

print_usage
exit 1
