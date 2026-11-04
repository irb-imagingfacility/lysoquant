# LysoQuant segmentation server with its own SSH (SFTP-only) server, for Fiji
# clients on other computers. Build context is the repository root: run
# ./setup.sh, then   docker compose up -d --build   (see README.md).
# PYTORCH_IMAGE and TORCH_INDEX are set by setup.sh (CUDA image or CPU image).
ARG PYTORCH_IMAGE=pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime
FROM ${PYTORCH_IMAGE}

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    LYSOQUANT_JOBS=/jobs \
    LYSOQUANT_WEIGHTS=/models/lysoquant-7-16.pt \
    LYSOQUANT_GPU=all \
    LYSOQUANT_SSH_PORT=2222

RUN apt-get update \
 && apt-get install -y --no-install-recommends openssh-server \
 && rm -rf /var/lib/apt/lists/* \
 && mkdir -p /run/sshd /etc/lysoquant

WORKDIR /opt/lysoquant
COPY server/pyproject.toml server/README.md ./
COPY server/lysoquant_server ./lysoquant_server
# torch comes with the CUDA base image; for the CPU image (python base) it is
# installed from TORCH_INDEX, the CPU-only wheels
ARG TORCH_INDEX=
RUN if [ -n "$TORCH_INDEX" ]; then pip install torch --index-url "$TORCH_INDEX"; fi \
 && pip install .

COPY model/lysoquant-7-16.pt /models/lysoquant-7-16.pt
COPY docker-ssh/sshd_config /etc/lysoquant/sshd_config
COPY docker-ssh/entrypoint.sh /usr/local/bin/lysoquant-entrypoint
RUN chmod 755 /usr/local/bin/lysoquant-entrypoint

EXPOSE 2222

# The segmentation server must be alive and sshd must listen. Looking in the
# socket table (state 0A = LISTEN) instead of connecting: a connection that
# closes before the SSH handshake is logged by sshd as an error
HEALTHCHECK --interval=30s --timeout=10s --start-period=120s --retries=3 \
    CMD ["bash", "-c", "lysoquant-health && grep -qiE \":$(printf '%04X' ${LYSOQUANT_SSH_PORT}) [0-9A-F:]+ 0A \" /proc/net/tcp /proc/net/tcp6"]

ENTRYPOINT ["lysoquant-entrypoint"]
