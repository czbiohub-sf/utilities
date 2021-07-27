FROM ubuntu:20.04

ENV PATH /opt/cellranger-arc-2.0.0:$PATH
ENV DEBIAN_FRONTEND=noninteractive

COPY cellranger-arc-2.0.0.tar.gz /opt

RUN ln -s /opt/cellranger-arc-2.0.0/bin/cellranger-arc /bin/cellranger-arc

RUN apt-get update && apt-get install -y git python3 python3-setuptools python3-pip

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN cd /opt && tar xvf cellranger-arc-2.0.0.tar.gz && rm cellranger-arc-2.0.0.tar.gz

RUN pip3 install --upgrade awscli

RUN cd /opt && \
    git clone https://github.com/czbiohub/utilities.git && \
    cd /opt/utilities && \
    git checkout neevor/cellrangerarc && \
    python setup.py install

WORKDIR /opt

CMD [ "/bin/bash"]