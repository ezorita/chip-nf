FROM ubuntu:14.04.5
RUN apt-get update -qq &&   \
    apt-get install -y      \
            build-essential \
            git             \
            wget            \
            gzip            \
            zlib1g-dev &&   \
    apt-get clean &&        \
    rm -rf /var/lib/apt/lists*
RUN cd / && \
    git clone http://github.com/lh3/bwa && \
    cd bwa && \
    git checkout 5961611c358e480110793bbf241523a3cfac049b && \
    make && \
    cp bwa /usr/bin
RUN cd / && \
    git clone http://github.com/nanakiksc/zerone && \
    cd zerone && \
    git checkout 449e8289244b38185f0ff37472d3ff55288d9615 && \
    make && \
    cp zerone /usr/bin
