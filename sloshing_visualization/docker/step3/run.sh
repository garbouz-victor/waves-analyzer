#!/usr/bin/env bash
set -euo pipefail
step3_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
step3_uid="$(id -u):$(id -g)"
# Docker Desktop's VirtioFS maps the host user to uid 0 INSIDE its VM.
# Native Linux Docker instead needs the host uid. No privileged container.
if [[ "$(docker context show)" == "desktop-linux" ]]; then
  step3_uid="0:0"
fi
exec docker run --rm --network=none --user "${STEP3_DOCKER_USER:-$step3_uid}" \
  --mount "type=bind,source=${step3_project_dir},target=/work" \
  waves-step3:0.10.0 "$@"
