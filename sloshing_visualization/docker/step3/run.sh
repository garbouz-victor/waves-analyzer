#!/usr/bin/env bash
set -euo pipefail
step3_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
step3_uid="$(id -u):$(id -g)"
# Docker Desktop's VirtioFS maps the host user to uid 0 INSIDE its VM.
# Native Linux Docker instead needs the host uid. No privileged container.
if [[ "$(docker context show)" == "desktop-linux" ]]; then
  step3_uid="0:0"
fi
step3_commit="$(git -C "$step3_project_dir" rev-parse HEAD)"
step3_image="$(docker image inspect --format '{{.Id}}' waves-step3:0.10.0)"
exec docker run --rm --network=none --user "${STEP3_DOCKER_USER:-$step3_uid}" \
  --env "STEP3_GIT_COMMIT=$step3_commit" --env "STEP3_DOCKER_IMAGE_DIGEST=$step3_image" \
  --env XDG_CACHE_HOME=/work/validation_results/step3a2/cache \
  --env MPLCONFIGDIR=/work/validation_results/step3a2/cache/matplotlib \
  --mount "type=bind,source=${step3_project_dir},target=/work" \
  waves-step3:0.10.0 "$@"
