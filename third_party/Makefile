


all: build

build:
	git submodule init && git submodule update && cd htslib/ && git checkout master && make && cd ../samtools/ && git checkout master && make && cd ../

clean:
	cd htslib/ && make clean && cd ../samtools/ && make clean

