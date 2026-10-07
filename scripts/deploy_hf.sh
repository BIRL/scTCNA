#!/usr/bin/env bash
# Deploy the current git HEAD to a Hugging Face Space (Docker SDK).
#
# Required environment:
#   HF_TOKEN  access token with write permission (https://huggingface.co/settings/tokens).
#             Read from the environment only: it is never written to disk, never stored in a
#             git remote, and only appears inline in the URL of the git commands that need it.
#   HF_SPACE  Space id, e.g. username/sctcna. Created with sdk=docker if it does not exist.
#
# The Space needs a README.md with Hugging Face YAML front matter, which the GitHub README
# should not carry. HEAD is therefore checked out into a temporary worktree, where
# scripts/hf_space_README.md is committed as README.md, and that commit is force-pushed to
# the Space's main branch. The current branch and working tree are left untouched.
#
# Prerequisites and manual fallbacks: docs/DEPLOY.md.
set -euo pipefail

: "${HF_TOKEN:?set HF_TOKEN to a Hugging Face write token}"
: "${HF_SPACE:?set HF_SPACE to the Space id, e.g. username/sctcna}"

case "$HF_SPACE" in
  */*) ;;
  *) echo "error: HF_SPACE must look like <owner>/<name>, got '$HF_SPACE'" >&2; exit 1 ;;
esac
OWNER="${HF_SPACE%%/*}"
NAME="${HF_SPACE#*/}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

SPACE_URL="https://huggingface.co/spaces/${HF_SPACE}"
# Used inline only; git redacts credentials from the URLs it prints.
AUTH_URL="https://${OWNER}:${HF_TOKEN}@huggingface.co/spaces/${HF_SPACE}"

# ---------------------------------------------------------------- preflight
# Everything below looks at HEAD, because that is what gets pushed.
if ! git rev-parse --verify -q HEAD >/dev/null; then
  echo "error: the repository has no commits yet" >&2
  exit 1
fi
HEAD_SHA="$(git rev-parse --short HEAD)"

for path in Dockerfile requirements-app.txt scripts/hf_space_README.md data/app; do
  if ! git cat-file -e "HEAD:$path" 2>/dev/null; then
    echo "error: $path is not in the committed tree (HEAD $HEAD_SHA); the Space builds from git." >&2
    [ "$path" = data/app ] && echo "Run 'make aggregate' and commit data/app first." >&2
    exit 1
  fi
done

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "note: uncommitted changes are ignored; only HEAD ($HEAD_SHA) is deployed" >&2
fi

# huggingface.co rejects plain git blobs over 10 MiB: such files must be git-lfs pointers.
# For LFS-tracked files the committed blob is the small pointer, so the size test passes.
limit=$((10 * 1024 * 1024))
problems=0
uses_lfs=0
while IFS= read -r f; do
  if git check-attr filter -- "$f" | grep -q 'filter: lfs$'; then
    uses_lfs=1
    continue
  fi
  size="$(git cat-file -s "HEAD:$f")"
  if [ "$size" -gt "$limit" ]; then
    echo "error: $f is $((size / 1024 / 1024)) MiB in HEAD and not tracked by git-lfs" >&2
    problems=1
  fi
done < <(git ls-tree -r --name-only HEAD)
if [ "$problems" -ne 0 ]; then
  echo "Track large files with git-lfs (e.g. git lfs track 'data/app/**/*.parquet'), commit, and retry." >&2
  exit 1
fi
if [ "$uses_lfs" -eq 1 ] && ! git lfs version >/dev/null 2>&1; then
  echo "error: the repository uses git-lfs but git-lfs is not installed (https://git-lfs.com)" >&2
  exit 1
fi

# Fail early on a bad token so a missing-Space error is not misdiagnosed.
# The header is piped to curl rather than passed as an argument.
if command -v curl >/dev/null 2>&1; then
  if ! printf 'Authorization: Bearer %s\n' "$HF_TOKEN" | curl -fsS -H @- https://huggingface.co/api/whoami-v2 >/dev/null; then
    echo "error: huggingface.co rejected HF_TOKEN (it needs write permission)" >&2
    exit 1
  fi
fi

# ------------------------------------------------- create the Space if missing
create_space() {
  if command -v hf >/dev/null 2>&1; then
    hf repos create "$HF_SPACE" --repo-type space --sdk docker
  elif command -v huggingface-cli >/dev/null 2>&1; then
    huggingface-cli repo create "$NAME" --organization "$OWNER" --type space --space_sdk docker -y
  else
    cat >&2 <<MSG
error: Space $HF_SPACE does not exist and no Hugging Face CLI is installed to create it.
Either install the CLI (pip install -U huggingface_hub) and re-run, or create it by hand at
https://huggingface.co/new-space (owner: $OWNER, name: $NAME, SDK: Docker), then re-run.
MSG
    return 1
  fi
}

if git ls-remote "$AUTH_URL" >/dev/null 2>&1; then
  echo "Space $HF_SPACE exists; updating it"
else
  echo "Space $HF_SPACE not found; creating it (sdk: docker)"
  if ! create_space; then
    echo "error: could not create the Space. Docker Spaces need a PRO account or a Team/Enterprise" >&2
    echo "organisation (see docs/DEPLOY.md); otherwise create it in the web UI and re-run." >&2
    exit 1
  fi
fi

# ------------------------------------------- push HEAD + Space README to main
TMP="$(mktemp -d "${TMPDIR:-/tmp}/sctcna-deploy.XXXXXX")"
cleanup() {
  cd "$ROOT" 2>/dev/null || true
  git worktree remove --force "$TMP" >/dev/null 2>&1 || rm -rf "$TMP"
}
trap cleanup EXIT

git worktree add --detach --quiet "$TMP" HEAD
cp scripts/hf_space_README.md "$TMP/README.md"

GIT_NAME="$(git config user.name || true)"
GIT_EMAIL="$(git config user.email || true)"
(
  cd "$TMP"
  git add README.md
  # The commit only ever lives on the Space, so it is not signed (no gpg dependency).
  if ! git diff --cached --quiet; then
    git -c user.name="${GIT_NAME:-sctcna deploy}" -c user.email="${GIT_EMAIL:-deploy@localhost}" \
      -c commit.gpgsign=false commit --quiet -m "Space README with Hugging Face front matter"
  fi
  # Each deploy is a fresh commit on top of HEAD, so the Space history is replaced, not extended.
  git push --force "$AUTH_URL" HEAD:main
)

APP_HOST="$(printf '%s-%s' "$OWNER" "$NAME" | tr 'A-Z_.' 'a-z--')"
cat <<MSG

Deployed $HEAD_SHA to $HF_SPACE
  Space page : $SPACE_URL
  Build logs : $SPACE_URL?logs=build
  App URL    : https://${APP_HOST}.hf.space   (live once the build finishes)
MSG
