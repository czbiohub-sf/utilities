FROM ubuntu:20.04

ENV PATH /opt/cellranger-arc-2.0.0:$PATH

COPY cellranger-arc-2.0.0.tar.gz /opt

RUN apt-get update && apt-get install -y git python3 python3-setuptools

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN cd /opt && tar xvf cellranger-arc-2.0.0.tar.gz && rm cellranger-arc-2.0.0.tar.gz

RUN cd /opt && \
    git clone https://github.com/czbiohub/utilities.git && \
    cd /opt/utilities && \
    git checkout neevor/cellrangerarc

WORKDIR /opt

CMD [ "/bin/bash"]