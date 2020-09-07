# Set the base image
FROM centos

#Load dependencies

RUN yum install -y \
git \
python2 \
wget \
epel-release \
python2-pip \
gcc \ 
python2-devel \
make \
zlib-devel \
gcc-c++ \
bzip2 \
bzip2-devel \
ncurses-devel \
xz-devel \
perl-Env \
java-devel

RUN yum install -y which

RUN dnf install -y redhat-rpm-config

#Get STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.4a.tar.gz
RUN tar -xzvf 2.7.4a.tar.gz
WORKDIR /STAR-2.7.4a/source
RUN make STAR
WORKDIR /
ENV PATH="/STAR-2.7.4a/bin/Linux_x86_64_static:${PATH}"
RUN STAR

#Get phaser:
RUN git clone https://github.com/secastel/phaser.git
RUN pip2 install Cython
RUN pip2 install scipy
RUN pip2 install pysam
WORKDIR phaser/phaser
RUN python2 setup.py build_ext --inplace
WORKDIR /

#Get Alleleseq:
RUN wget http://alleleseq.gersteinlab.org/vcf2diploid_v0.2.6a.zip && unzip vcf2diploid_v0.2.6a.zip
WORKDIR vcf2diploid_v0.2.6a
RUN make
WORKDIR /

#Get Samtools:
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN bzip2 -d samtools-1.10.tar.bz2
RUN tar -xvf samtools-1.10.tar
RUN echo $(ls)
WORKDIR samtools-1.10
RUN ./configure --prefix=/usr
RUN make
RUN make install
WORKDIR /usr/bin
RUN echo $(ls)
WORKDIR /

#Liftover:
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
RUN echo $(ls)
RUN chmod +x liftOver && mv liftOver /usr/bin

#RSEM
RUN git clone https://github.com/deweylab/RSEM.git
WORKDIR RSEM
RUN make
WORKDIR /
ENV PATH="/RSEM:${PATH}"

#Get HTSLIB: 
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 
RUN bzip2 -d htslib-1.10.2.tar.bz2 
RUN tar -xvf htslib-1.10.2.tar 
RUN echo $(ls) 
WORKDIR htslib-1.10.2 
RUN ./configure --prefix=/usr 
RUN make 
RUN make install 
WORKDIR /usr/bin 
RUN echo $(ls) 
WORKDIR /

