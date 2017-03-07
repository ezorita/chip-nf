FROM ubuntu:14.04.5
RUN apt-get update -qq &&   \
    apt-get install -y      \
            build-essential \
            git             \
            wget            \
            gzip            \
            zlib1g-dev      \
            libxml2-dev     \
            libmagic-dev    \
            libhdf5-dev     \
            libfuse-dev &&  \
    apt-get clean &&        \
    rm -rf /var/lib/apt/lists*
RUN cd / && \
    git clone http://github.com/lh3/bwa && \
    cd bwa && \
    git checkout 5961611c358e480110793bbf241523a3cfac049b && \
    make && \
    cp bwa /usr/local/bin
RUN cd / && \
    git clone http://github.com/nanakiksc/zerone && \
    cd zerone && \
    git checkout 449e8289244b38185f0ff37472d3ff55288d9615 && \
    make && \
    cp zerone /usr/local/bin
RUN cd / && \
    git clone http://github.com/ezorita/bioscripts && \
    cd bioscripts && \
    git checkout 3576cfd3f9a77863bd59899eaf2a2ccf5edfd148 && \
    make colortobase && \
    cp colortobase /usr/local/bin
RUN cd / && \
    mkdir ncbi && \
    cd ncbi && \
    git clone http://github.com/ncbi/ngs && \
    git clone http://github.com/ncbi/ncbi-vdb && \
    git clone http://github.com/ncbi/sra-tools && \
    ./ngs/configure && make -C ngs/ngs-sdk install && \
    ./ncbi-vdb/configure && make -C ncbi-vdb install && \
    ./sra-tools/configure && make -C sra-tools install && \
    mv /usr/local/ncbi/sra-tools/bin/* /usr/local/bin/
