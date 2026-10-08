# Deploying the scTCNA web app

The app is a FastAPI service that serves the precomputed Parquet tables in
`data/app` (built by `python -m sctcna aggregate`). Nothing is recomputed at
request time, so it runs in a small container: Python 3.12, FastAPI, uvicorn,
pandas, pyarrow and numpy. The pipeline stack (scanpy, anndata, Scrublet,
SciPy, statsmodels) is not installed in the image.

| File | Purpose |
| --- | --- |
| `Dockerfile` | `python:3.12-slim`, unprivileged user, installs `requirements-app.txt`, copies `sctcna/`, `app/`, `configs/`, `data/app/`; serves on `$PORT` (default 7860) |
| `requirements-app.txt` | runtime pins (a subset of `requirements.txt`) |
| `render.yaml` | Render Blueprint: one free Docker web service with a health check on `/api/health` |
| `.dockerignore` | keeps raw data, per-dataset outputs, tests and caches out of the build context |

## Host: Render free tier

Render's free web services cost nothing and need no payment card. Limits that
matter here: 512 MB RAM, a shared CPU, and the service sleeps after 15 minutes
without traffic (the next request wakes it in roughly a minute).

The app was measured against the full ten-cancer tables before choosing this
host; the numbers are in the "Memory" section below. Hugging Face Docker Spaces
were ruled out because creating one now requires a paid plan.

## Steps

1. **Make sure `data/app` is committed.** The image copies it, and Render builds
   from the repository. It is committed in this repository; after re-running
   the pipeline, run `python -m sctcna aggregate` and commit `data/app` again.
2. **Push the repository to GitHub** (public or private):

   ```bash
   git remote add origin git@github.com:<your-user>/sctcna.git
   git push -u origin main
   ```

3. **Create the service from the Blueprint.**
   - Sign in at <https://dashboard.render.com> (GitHub login is fine; the free
     plan does not ask for a card).
   - Click **New** > **Blueprint**, connect the GitHub account if asked, and
     pick the `sctcna` repository.
   - Render reads `render.yaml` and shows one web service, `sctcna`, plan
     *Free*. Click **Apply**.
4. **Wait for the first build** (about 5 minutes: the image installs five
   wheels). The Logs tab ends with `Uvicorn running on http://0.0.0.0:...`
   and the service turns *Live*.
5. **Open the URL** shown at the top of the service page
   (`https://sctcna-XXXX.onrender.com`) and check `/api/health`; it should
   report `"status":"ok"` and `"n_datasets":10`.

Every later `git push` to `main` redeploys automatically (`autoDeploy: true`).

Without the Blueprint, the same service can be created by hand: **New** >
**Web Service** > pick the repository > *Language*: Docker > *Instance type*:
Free > *Advanced*: health check path `/api/health`, environment variable
`SCTCNA_DATA_DIR=/app/data/app` > **Deploy**.

## Memory

Measured locally with `/usr/bin/time -l` on a uvicorn process serving the
committed `data/app` (ten cancer types) after requesting every page and the
heaviest endpoints (all comparisons, the per-cell violin data of every cancer
type, a 200-gene heatmap). See `docs/memory_check.txt` for the exact numbers
of the last check; the release was checked to stay below 512 MB with margin.

Three things keep it small: the app tables are stored with categorical strings
and float32 display values, the per-patient tables are loaded on demand with at
most three in memory, and the app does not import SciPy or statsmodels (the
paired DEPTH tests are precomputed by `sctcna aggregate`).

## Run the container locally

```bash
docker build -t sctcna .
docker run --rm -p 7860:7860 sctcna
```

Then open <http://localhost:7860>. To try the image with other tables, mount
them over the baked-in directory:

```bash
docker run --rm -p 7860:7860 -v "$PWD/data/app_syn:/app/data/app:ro" sctcna
```

| Variable | Default | Meaning |
| --- | --- | --- |
| `SCTCNA_DATA_DIR` | `/app/data/app` | directory with the Parquet tables and `build.json` |
| `PORT` | `7860` | port uvicorn binds to; Render sets it automatically |

## Troubleshooting

- **Build fails at `COPY data/app/`**: the tables are not in the repository.
  Run `python -m sctcna aggregate` and commit `data/app`.
- **Service restarts with "out of memory"**: re-check with the procedure in
  `docs/memory_check.txt`; reduce `MAX_PATIENT_TABLES` in `app/store.py` to 1.
- **First request after a pause is slow**: the free service was asleep; it
  wakes on demand.
