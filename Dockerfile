FROM arkadiuszgladki/gdr_shiny:0.08


# GitHub token for downloading private dependencies
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY gDRtestData/ /tmp/gDRtestData/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRtestData')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRtestData/
