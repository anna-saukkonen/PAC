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
java-devel

RUN dnf install -y redhat-rpm-config

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

