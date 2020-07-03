# Set the base image to debian jessie
FROM debian:jessie

# File Author / Maintainer
MAINTAINER Anna Saukkonen <anna.saukkonen@gmail.com>



RUN apt-get update && apt-get install --yes --no-install-recommends \
	ca-certificates \
    wget \
    curl \
    locales \
    vim-tiny \
    git \
    cmake \
    build-essential \
    gcc-multilib \
    perl \
    python

RUN wget https://github.com/secastel/phaser.git \
	&& cd phaser/phaser/ \
	&& python setup.py build_ext –inplace  

RUN wget http://alleleseq.gersteinlab.org/vcf2diploid_v0.2.6a.zip \
	&& make	
