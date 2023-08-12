# Use a base image with Ubuntu
FROM ubuntu:latest

# Set the working directory
WORKDIR /app

# Install required packages
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        g++ \
        wget \
        automake \
        libtool \
        git
# Copy the Spot library tar.gz file into the Docker image
COPY spot-2.11.6.tar.gz /app

# Extract and install Spot
RUN tar -xf spot-2.11.6.tar.gz && \
    cd spot-* && \
    ./configure --disable-python && \
    make && \
    make install

# Install CUDD
RUN git clone https://github.com/KavrakiLab/cudd.git && \
    cd cudd && \
    autoreconf -i && \
    ./configure --enable-silent-rules --enable-obj --enable-dddmp --prefix=/usr/local && \
    make && \
    make install

# Install MONA
RUN apt-get update && apt-get install -y flex && \
    wget https://www.brics.dk/mona/download/mona-1.4-17.tar.gz && \
    tar -xf mona-1.4-17.tar.gz && \
    cd mona-* && \
    ./configure && \
    make && \
    make install-strip

# Set environment variables
ENV MONA_BIN=/usr/local/bin/mona
ENV SPOT_BIN=/usr/local/bin/ltl2tgba
ENV CUDD_DIR=/usr/local
ENV LD_LIBRARY_PATH /usr/local/lib

# Copy your source code
COPY . /app

# Include Spot include path
ENV CPLUS_INCLUDE_PATH=/usr/local/include/spot

# Compile your program
RUN g++ -std=c++17 lisa.cc dfwavar.cc dfwa.cc spotutil.cc ltlf2fol.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc -o lisa -lspot -lbddx -lcudd -O3

# Set the default command to run your program
CMD ["./lisa"]
