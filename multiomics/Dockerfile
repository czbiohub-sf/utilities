FROM ubuntu:20.04

ENV PATH /opt/cellranger-arc-2.0.0:/opt/cellranger-6.1.0:$PATH
ENV DEBIAN_FRONTEND=noninteractive

COPY cellranger-arc-2.0.0.tar.gz /opt
COPY cellranger-6.1.0.tar.gz /opt

RUN ln -s /opt/cellranger-arc-2.0.0/bin/cellranger-arc /bin/cellranger-arc
RUN ln -s /opt/cellranger-6.1.0/bin/cellranger /bin/cellranger

RUN apt-get update && apt-get install -y git python3 python3-setuptools python3-pip

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN cd /opt && \
    tar xvf cellranger-arc-2.0.0.tar.gz && \
    rm cellranger-arc-2.0.0.tar.gz && \
    tar xvf cellranger-6.1.0.tar.gz && \
    rm cellranger-6.1.0.tar.gz

RUN pip3 install --upgrade awscli \
                           pandas==1.3.1 \
			   "numpy<1.21,>=1.17" \
                           scanpy==1.8.1

RUN cd /opt && \
    git clone https://github.com/czbiohub/utilities.git && \
    cd /opt/utilities && \
    git checkout neevor/cellrangerarc && \
    python setup.py install

WORKDIR /opt

CMD [ "/bin/bash"]
