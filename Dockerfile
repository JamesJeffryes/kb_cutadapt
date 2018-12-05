FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer

# Here we install a python coverage tool
RUN pip install coverage

# install cutadapt
RUN pip install cutadapt==1.18

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod 777 /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
