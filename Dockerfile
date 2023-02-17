ARG BASE_IMAGE=marcinkam/gdrshiny:0.11
FROM ${BASE_IMAGE}

# GitHub token for downloading private dependencies
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN echo "$GITHUB_TOKEN" >> /mnt/vol/.github_access_token.txt
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY gDRtestData/ /tmp/gDRtestData/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRtestData')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRtestData/
