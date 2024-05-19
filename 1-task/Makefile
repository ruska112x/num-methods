project = numeric1

all: build exec

run: reload exec

build: mkdir-build cmake-load cmake-build

reload: cmake-load cmake-build

mkdir-build:
	[ -d ./build ] | mkdir -p build

cmake-load:
	cd build;cmake ..

cmake-build:
	cd build;cmake --build . --target $(project)

clean:
	rm -rf build

exec:
	./build/$(project)
