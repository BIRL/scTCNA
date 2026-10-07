# Deploying the scTCNA web app

The app is a FastAPI service that reads the precomputed Parquet tables in
`data/app` (built by `make pipeline aggregate`). Nothing is recomputed at
request time, so the deployment is a small container: Python 3.12, FastAPI,
uvicorn, pandas, pyarrow, numpy, scipy and statsmodels. The single-cell
pipeline stack (scanpy, anndata, scrublet, scikit-learn) is not installed in
the image.

Files involved:

| File | Purpose |
| --- | --- |
| `Dockerfile` | `python:3.12-slim`, non-root user `user` (uid 1000), installs `requirements-app.txt`, copies `sctcna/`, `app/`, `configs/`, `data/app/`, `pyproject.toml`; serves on port 7860 |
| `.dockerignore` | keeps raw data, per-dataset outputs, synthetic tables, tests and caches out of the build context |
| `requirements-app.txt` | runtime pins (must match `requirements.txt`) |
| `scripts/deploy_hf.sh` | pushes git HEAD to a Hugging Face Space, creating it if needed |
| `scripts/hf_space_README.md` | README with the Spaces YAML front matter; copied to `README.md` inside the Space only |
| `Makefile` | `docker-build`, `docker-run`, `deploy-hf` and the development targets |

## Where it runs and why

Target: **Hugging Face Spaces, Docker SDK, CPU Basic hardware** (2 vCPU,
16 GB RAM, 50 GB non-persistent disk, no hourly charge, persistent URL
`https://<owner>-sctcna.hf.space`, automatic rebuild on every push).

Why this over the other no-cost hosts:

- **Render free web services** give 512 MB RAM and spin down after 15 minutes
  of inactivity. pandas + pyarrow with the atlas tables resident in memory is
  tight at 512 MB and every cold start reloads them; 16 GB on Spaces leaves
  room to spare.
- **Railway** only offers a time-limited trial credit, after which the service
  is paid.
- Spaces keep the code and the deployment together in one git repository,
  build the Dockerfile themselves (no registry to push to) and need no
  payment method for the CPU Basic hardware itself.

Caveat, checked against the Hugging Face documentation on 2026-10-08
(<https://huggingface.co/docs/hub/spaces-overview>): since mid-2026 **creating
a Gradio or Docker Space requires a paid plan** (PRO for personal accounts,
Team or Enterprise for organisations); only static Spaces are free to create.
The CPU Basic hardware is still listed at no hourly cost, so the running cost
is the PRO subscription, not compute. If that is not acceptable, the same
image runs unchanged on any container host that injects a `PORT` variable
(Render, Railway, Fly.io, Google Cloud Run, a VM with Docker); the Dockerfile
honours `PORT` for that reason. The choice of host is the owner's call; this
document describes the Spaces route.

Free-hardware Spaces go to sleep after 48 hours without traffic and wake on
the next request (a cold start takes about as long as the container start,
a few seconds, because the tables are small).

## Prerequisites

1. **`data/app` committed to git.** The Space builds from the repository, and
   `COPY data/app/ ./data/app/` fails if the directory is missing. Build it
   with `make pipeline aggregate` (or `make aggregate` if `data/processed`
   already exists) and commit the result. `.gitignore` already keeps
   `data/raw` and `data/processed` out and explicitly allows `data/app`.
2. **Files over 10 MiB must be git-lfs pointers.** huggingface.co rejects
   pushes containing plain git blobs larger than 10 MiB ("Your push was
   rejected because it contains files larger than 10 MiB"). `cell_scores.parquet`
   (one row per cell) and the per-patient tables will exceed that on real
   datasets. Either split them or track them with LFS:

   ```bash
   git lfs install
   git lfs track 'data/app/**/*.parquet'
   git add .gitattributes data/app
   git commit -m "Track app tables with git-lfs"
   ```

   The deploy script checks this and refuses to push otherwise. LFS objects
   also count towards GitHub's Git LFS storage and bandwidth quota, which on
   free plans is ample for a few hundred megabytes of tables.
3. A Hugging Face account on a plan that can create Docker Spaces (see the
   caveat above) and an access token with **write** permission from
   <https://huggingface.co/settings/tokens>.
4. Optional: the Hugging Face CLI (`pip install -U huggingface_hub` gives
   `hf`; older versions ship `huggingface-cli`). The script uses it only to
   create a missing Space; an existing Space needs git alone.

## Deploy with the script

```bash
export HF_TOKEN=<write token>          # from the environment only; never put it in a file
export HF_SPACE=<owner>/sctcna         # user or organisation namespace
make deploy-hf                          # or: bash scripts/deploy_hf.sh
```

What it does, in order:

1. Checks that `Dockerfile`, `requirements-app.txt`,
   `scripts/hf_space_README.md` and `data/app` are in the committed HEAD
   tree, that no committed file over 10 MiB lacks an LFS pointer, and (when
   `curl` is available) that the token is accepted by huggingface.co.
2. Looks for the Space with `git ls-remote`. If it is missing, creates it
   with `hf repos create <space> --repo-type space --sdk docker` (or the
   `huggingface-cli` equivalent); without a CLI it prints the manual steps
   and stops.
3. Checks HEAD out into a temporary detached worktree, copies
   `scripts/hf_space_README.md` over `README.md` there and commits it, so the
   Space gets the YAML front matter while the GitHub README stays as it is.
4. Force-pushes that commit to the Space's `main` branch with the token in
   the push URL of that single command. No remote is added and nothing is
   written to disk; the worktree is removed on exit.
5. Prints the Space page, the build-log URL and the app URL.

Only the committed HEAD is deployed; uncommitted changes are ignored (the
script says so). Re-run the same command to redeploy after new commits. The
Space rebuilds automatically after each push; a build takes a few minutes
(pip install of the pinned wheels dominates).

Visibility defaults to public. Change it under the Space's **Settings**
tab if needed; a private Space still serves the app to its collaborators.

## Manual fallbacks

### A. Push with git by hand

```bash
# 1. Create the Space once: https://huggingface.co/new-space
#    Owner: <owner>, Space name: sctcna, SDK: Docker, hardware: CPU basic.
# 2. From the repository root, on the commit you want to deploy:
git worktree add --detach /tmp/sctcna-space HEAD
cp scripts/hf_space_README.md /tmp/sctcna-space/README.md
git -C /tmp/sctcna-space add README.md
git -C /tmp/sctcna-space commit -m "Space README with Hugging Face front matter"
git -C /tmp/sctcna-space push --force "https://<owner>:$HF_TOKEN@huggingface.co/spaces/<owner>/sctcna" HEAD:main
git worktree remove --force /tmp/sctcna-space
```

Git prompts for credentials if the token is left out of the URL: use your
username and the token as the password. Do not run `git remote add` with the
token in the URL, or it ends up in `.git/config`.

### B. Upload through the web UI (no git, no CLI)

1. Open <https://huggingface.co/new-space>. Owner: your account or
   organisation; Space name: `sctcna`; License: MIT; SDK: **Docker** (blank
   template); hardware: CPU basic; visibility as you prefer. Click
   **Create Space**.
2. In the new Space open the **Files** tab, then **Contribute > Upload files**.
3. Drag in, from the repository root: `Dockerfile`, `requirements-app.txt`,
   `pyproject.toml`, and the folders `sctcna/`, `app/`, `configs/` and
   `data/app/` (the uploader keeps folder structure when you drop a folder;
   files over 10 MiB are stored with LFS automatically). Add a commit message
   and click **Commit changes to main**.
4. Open `README.md` in the Files tab, click the edit (pencil) icon, replace
   the whole content with `scripts/hf_space_README.md` and commit. The
   front matter must contain `sdk: docker` and `app_port: 7860`.
5. The **App** tab shows the build log; the app is live when the status badge
   turns to *Running*.

To update later, upload the changed files the same way (uploading a file with
an existing path replaces it).

## Run locally with Docker

```bash
make docker-build                 # docker build -t sctcna:latest .
make docker-run                   # docker run --rm -p 7860:7860 sctcna:latest
open http://localhost:7860        # or curl -s localhost:7860/api/health
```

Both need `data/app` to exist in the working tree. To try the image with the
synthetic tables instead, mount them over the baked-in directory:

```bash
docker run --rm -p 7860:7860 -v "$PWD/data/app_syn:/app/data/app:ro" sctcna:latest
```

Runtime configuration (all optional):

| Variable | Default | Meaning |
| --- | --- | --- |
| `SCTCNA_DATA_DIR` | `/app/data/app` | directory with the Parquet tables and `build.json` |
| `PORT` | `7860` | port uvicorn binds to; Spaces expect 7860 (`app_port` in the README) |

The image runs one uvicorn worker, which is enough for read-only Parquet
queries on 2 vCPUs and keeps one copy of the tables in memory. The container
has a `HEALTHCHECK` on `/api/health` (used by `docker run`/compose; Spaces
ignore it and rely on the app answering on `app_port`).

Image size: `python:3.12-slim` (about 130 MB) plus the pinned wheels (pyarrow
and scipy are the large ones) comes to roughly 550-650 MB uncompressed before
`data/app`; installing the full `requirements.txt` would roughly double that,
which is why `requirements-app.txt` exists.

## Troubleshooting

- **Build fails at `COPY data/app/`**: the directory is not in the Space
  repository. Commit `data/app` (prerequisite 1) and push again.
- **Push rejected, "files larger than 10 MiB"**: see prerequisite 2.
- **`hf repos create` or the web UI refuses a Docker Space**: the account is
  on a free plan (see the caveat above). Upgrade, use an organisation on a
  Team plan, or deploy the same image elsewhere.
- **Space shows *Runtime error* with the container exiting immediately**:
  open the build and container logs (`?logs=container` on the Space URL).
  The usual cause is a missing table: the app expects every file listed in
  the data contract under `SCTCNA_DATA_DIR`.
- **App is reachable but slow on first request after idle**: the Space was
  asleep; it wakes on demand on free hardware.
- **Updating pins**: change `requirements.txt` and `requirements-app.txt`
  together, rebuild locally with `make docker-build`, then deploy.
