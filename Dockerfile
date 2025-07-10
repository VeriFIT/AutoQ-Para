FROM fedora:latest
RUN dnf install -y cmake make g++ gmp-devel
RUN mkdir -p /workspace/build
COPY ./src ./workspace
WORKDIR /workspace/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j6 experiments
ENTRYPOINT ./experiments
