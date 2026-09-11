# Проверенная база и что переиспользовать

Проверено чтением GitHub 2026-09-10:
- feature/pinned-wetting-history -> ac418484d48dff5f8104d2177b5c2b1f703a232a;
- master -> 7f0c45f192520b218838d535b0182facf43133d3.
Это наблюдение, не приказ переставлять ветки. Перепроверить локальные changes.

`src/sloshing/pinned_wetting/asymmetric.py`: правильные left Navier/right+bottom
no-slip, trace [1,0]; это LINEAR fixed-domain FEM, не moving-domain solver.
`asymmetric_run.py`: identity/provenance и bridge выполнения.
`asymmetric_verify.py`: independent left-only BC assembly, geometry witnesses.
`asymmetric_mission.py`, `asymmetric_package.py`: pilot, отчёт и явно failed release.
`coating.py`: постоянное связное покрытие и immutable markers.
`mission.py`: job lock, бюджет/heartbeat, CLI. Не создавать параллельный ledger.
`extension_output.py`: существующий renderer, нужно адаптировать на Gamma, а не
просто скормить немонотонную границу функции z=eta(x).
`extension_verify.py`: полезный пример независимых stage/work проверок; его
линейные формулы не являются доказательством nonlinear correctness.
`multiphase/full_phase_rate.py`: ограничен matched density/zero force; не применять
к gravity target отключением require_scope. Не обязан завершать STEP3A ради новой задачи.

Текущий corrected pilot: правильная P2, но slope max 9.44/15.70 и global surface
mesh discrepancy 18.24%. Это основание менять приближение/представление, не BC.
Из старых reports не следует, что sigma=0 nonlinear проблема неразрешима, и
не следует, что одна замена координат обязательно вылечит все особенности.

Нативные HDF5 и MP4 в основном local ignored. Работай в существующем worktree
с этими данными; отсутствующие файлы не объявлять прочитанными/проверенными.
