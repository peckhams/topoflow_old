FROM centos:7

LABEL maintainer="Scott D. Peckham <Scott.Peckham@colorado.edu>" contributor="Rajiv Mayani <mayani@isi.edu>"

RUN yum -y update && yum -y install epel-release && yum -y install python-setuptools python-pip numpy scipy
# RUN yum -y update && yum -y install epel-release && yum -y install python36

RUN mkdir /srv/topoflow

ADD . /srv/topoflow

RUN cd /srv/topoflow && pip install .
