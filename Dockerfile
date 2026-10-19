# Runtime image for the scTCNA web app: FastAPI serving the precomputed tables
# in data/app. The single-cell pipeline (scanpy and friends) is not installed;
# tables are built offline with `make pipeline aggregate` and committed.
# Deployed on Render's free tier (render.yaml); runs on any container host.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Run as an unprivileged user; creating it first lets every COPY below be owned
# by it without a separate chown layer.
RUN useradd -m -u 1000 user
WORKDIR /app

# Dependencies first so code changes do not invalidate the pip layer.
COPY requirements-app.txt ./
RUN pip install --no-cache-dir -r requirements-app.txt

COPY --chown=user:user pyproject.toml ./
COPY --chown=user:user sctcna/ ./sctcna/
COPY --chown=user:user app/ ./app/
COPY --chown=user:user configs/ ./configs/
# Precomputed app tables. data/app must exist in the build context (committed
# to git); see docs/DEPLOY.md.
COPY --chown=user:user data/app/ ./data/app/

# The system allocator returns freed Arrow buffers to the OS, which keeps resident memory low
ENV SCTCNA_DATA_DIR=/app/data/app \
    ARROW_DEFAULT_MEMORY_POOL=system \
    PORT=7860 \
    PYTHONPATH=/app \
    HOME=/home/user

USER user
EXPOSE 7860

# Render probes /api/health itself (render.yaml); `docker run` uses this check.
# The slim image has no curl, so the probe uses the standard library.
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD python -c "import os, urllib.request; urllib.request.urlopen('http://127.0.0.1:%s/api/health' % os.environ.get('PORT', '7860'), timeout=4)"

# Shell form so ${PORT} is honoured on hosts that inject their own port (Render
# sets PORT); exec makes uvicorn PID 1 so it receives SIGTERM.
CMD ["sh", "-c", "exec uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-7860} --workers 1"]
