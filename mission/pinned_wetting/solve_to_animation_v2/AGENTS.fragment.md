<!-- pw2-fixed-conditions:start -->
## Current mission: PW2 fixed conditions to final animation
Read `mission/pinned_wetting/solve_to_animation_v2/MASTER_PROMPT.md` and CONTRACT.json.
For this new mode, preserve LEFT Navier / RIGHT+BOTTOM no-slip / sigma=0.
The new prompt explicitly permits a separately labelled full nonlinear free-boundary
solve. Never ask the user to choose the predicted interface shape. Keep old linear
and symmetric results historical. Reuse existing CLI/locks/budget/checkpoints.
Continue through milestones; no fsync detours. Export guards are not PDE validation.
<!-- pw2-fixed-conditions:end -->
