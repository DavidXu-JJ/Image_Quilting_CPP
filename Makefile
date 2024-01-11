
.PHONY:clean
clean:
	rm -rf build

.PHONY:cmake
cmake:
	cmake -S . -B build

.PHONY:make
make:
	cmake --build build -j12

.PHONY:all
all:
	make clean
	make cmake
	make make

.PHONY: run
run:
	./build/image_quilting $(ARGS)
