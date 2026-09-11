# Измеримый финиш

## Название результата

`COMPLETE_FOR_DECLARED_FREE_BOUNDARY_MODEL` допустим только для target из
CONTRACT.json. Рядом явно указать full nonlinear NS и отличие от прежнего
линейного reference. Не переименовывать прошлые FAIL/COMPLETE.

## Обязательный результат

Новые настоящие данные 0–5 физических секунд, 11 срезов через 0.5 с,
`final_animation.mp4`, `wall_detail.mp4`, `index.html`. Поверхность в обоих
роликах из accepted native geometry. Полный просмотр пользователем не требует
HDF5 viewer: index содержит видео, ключевые кадры, таблицу результатов и отчёт.

Output root: `sloshing_visualization/output/pinned_wetting/free_boundary_fixed_conditions/`.
Минимальные сопутствующие файлы:

    resolved_case.json
    model_notes.md
    interface_samples.csv
    contacts.csv
    retained_intervals.csv
    marker_catalog.csv
    energy_history.csv
    refinement.json
    native_manifest.json
    native_verification.json
    review.md
    key_frames.png
    FINAL_REPORT.md
    manifest.json
    verification.json

Нативная геометрия может быть parametric/XDMF/HDF5/другой формой; это записано
в native_manifest. Не преобразовывать многозначную кривую в один eta, теряя ветвь.
Полный интерфейс сохранён на каждом accepted состоянии или восстанавливается
без неоднозначности; contact/history на каждом accepted шаге, а не только в видео.

## Обязательные инварианты

- Полный физический контракт совпадает с CONTRACT.json; новое sigma/right slip
  или изменение alpha/nu для PASS запрещено.
- Начальная поверхность — заданная прямая; v0=0. P1/P2 фиксированы как маркеры.
- Материальный RIGHT контакт сохраняет (+a,z2), velocity trace RIGHT/BOTTOM=0,
  нормальная скорость LEFT=0; tangential Navier LEFT проверен операторно.
- Нет правого Robin term даже в до-редуцированной матрице. Нет endpoint clamp
  в renderer. Поворот normal/касательного базиса проверен отдельно.
- H соответствует history реального contact detection, не просто монотонен;
  W=[-d,H] хранится без dry holes в checkpoint и переживает restart.
- Rejected steps не изменяют accepted state, H, маркеры или simulation time.
- Масса сохраняется; полные stress/kinematics/weak incompressibility проверены
  независимыми residuals, не только production summary.

Для exact Dirichlet/endpoint целью остаётся 1e-12 в SI. Для другого способа
наложения ограничений необходима обоснованная методическая оценка до production,
а не post-hoc увеличение допуска. Volume relative ≤1e-6; continuous-energy budget
relative ≤0.05 на E0, с отдельно измеренной численной диссипацией/ошибкой, а не
подогнанным закрывающим членом. Слишком слабый energy gate не заменяет PDE residual.

## Точность, которую можно повышать

Отдельные сравнения dt и h при одних BC/физике до 5 с. В двух уровнях formal order
не доказан. Цель: максимум геометрического расстояния полной поверхности /
|z2-z1|≤0.05; left contact и H differences/|z2-z1|≤0.05. Отдельно peak heights/times
и energy observables; near-threshold record identity может быть неопределённа.
Если нужно, одно дополнительное направленное локальное refinement в той же миссии.

Точный молекулярный film thickness не нужен. Но исключить ещё неразрешённый
физически важный участок границы из сравнения и назвать всю анимацию точной нельзя.
Большой slope НЕ automatic FAIL полной геометрии: geometry validity, Jacobian,
gap resolution, conserved transfer и refinement определяют достоверность.

## Выходы verifier

`contract_guard.py` и `check_exports.py` из пакета — только дешёвые внешние проверки.
Они не дают CFD PASS. Итоговый `verification.json` обязан разделять:

    contract_passed
    independent_native_equations_passed
    boundary_conditions_passed
    geometry_resolution_passed
    mass_energy_passed
    spatial_refinement_passed
    temporal_refinement_passed
    retained_coating_passed
    video_decode_passed
    independent_visual_review_passed

Каждый true ссылается на реальный native/расчётный evidence, digest и версию
проверяющего кода. Отдельный контекст review полезен, но не заменяет независимую
пересборку уравнений. Не называть самоподпись JSON независимой экспертизой.

## Негативные тесты

Минимум: right freed/right Robin добавлен, перепутаны стороны, изменён b,
sigma>0, fake moving-right-contact при неподвижной P2 в каталоге, ошибочная
normal, неверная ALE convection, нарушенный GCL, утечка при remesh, уменьшение H,
разрыв W, потеря истории на restart, snapshot другой physics, fake 5s из 1s,
линейный reference в nonlinear package, данные не от solver, потерянная тонкая
ветвь при экспорте, post-commit compatible resume. Отсутствие P4 — не ошибка.

## Видео и usability

0–5 с, физическое время на каждом кадре; 1920×1080,25fps и заявленное замедление,
если эти settings совместимы с доступным renderer. Плотность render samples
не определяет detector sampling. Оси/масштабы фиксированы и подписаны. Полный
сосуд и увеличенная пристеночная зона. Условная полоска смачивания визуально
отличается от рассчитанной разрешённой плёнки; не дорисовывать правую геометрию.
Маркеры не двигаются и не возникают до подтверждённого контакта/пика.

Видео декодировать полностью, проверить кадры t=0, первого left peak, t=1,2.5,5
и момента максимальной деформации справа. Проверить совпадение с accepted native.

## Честное незавершение

Статусы RESOURCE_LIMIT, UNRESOLVED_REGION, MODEL_CLOSURE_EVIDENCE не являются
COMPLETE. Они сопровождаются пригодным preview, реальным пределом времени и
диагностикой. Не требуют от пользователя заранее выбрать вид поверхности.
Гарантировать 5 с без расчёта нельзя; объявить же общий FAIL по старому slope
или cost forecast и остановиться без работы над методом тоже нельзя.
