# Карта существующего проекта — прочитать до выбора solver

Проверенный срез master: `7f0c45f192520b218838d535b0182facf43133d3`.
Проверка сделана 2026-09-10. При более новом HEAD выполнить короткий gap-audit,
не reset. Рекомендуемая новая ветка: `feature/pinned-wetting-history`.

## Переиспользовать

`sloshing_visualization/src/sloshing/config.py`: параметры в SI, a/d/nu/alpha,
mesh presets, integrator selection. Defaults nu=0.01, alpha=2° — это кодовые
defaults, не установленные пользователем свойства материала.

`fem_spaces.py`, `problems.py`, `time_integrator.py`, `solver.py`: линейная FEM
задача, начальная прямая, no-slip и численное интегрирование. Отличный regression
reference для исходного малого bulk sloshing; не новый wetting solver.

`src/sloshing/visualization/`, `scripts/render_step2.py`: video rendering,
фиксированные шкалы, no-slip diagrams, offline explorer. Проверить ожидания
формата dataset; не переименовывать старый HDF5 в target run.

`src/sloshing/multiphase/`: finite-element/phase-rate, energy, interface,
checkpoint и диагностика — полезные компоненты. Подключать только то, что
соответствует нужной постановке; не тащить весь historical pipeline как dependency.

## Ограничения, которые агент не имеет права скрыть

### Линейный fixed-domain endpoint

В `fem_spaces.py` wall velocity DOFs удалены из свободных. R отображает
вертикальную скорость на surface trace, включая endpoints; их columns уже
нулевые. В `time_integrator.py` eta обновляется через R v.
Значит eta endpoints сохраняют исходную высоту. Running maximum от eta[0]/eta[-1]
даст исходные значения. Near-wall extrapolation не создаёт физическое wall contact.

### Поздний full phase-rate CHNS

В `multiphase/full_phase_rate.py` `require_scope()` явно требует:
matched density, g=0, a_x=0, BE, P2, quadrature=12.
Существующий PASS нельзя перенести на gravity/unequal density удалением проверки.
Нужно ограниченное обоснование и тест действительно используемого нового пути.

### Film helper

`multiphase/film.py` содержит классификацию разрешимости и формулу Nusselt flux;
это не готовый coupled solver всей плёнки. Для этой миссии точная толщина не нужна.

## Ветки

`validation/step2-1-no-slip` — исторический visualization checkpoint.
`linear-streamfunction-semi-infinite` — отдельный linear formulation/design;
не замена готовому движению контакта. Использовать их как источники отдельных
проверок при необходимости, не переключать mission между ветками без нужды.

## Большие данные

MP4, некоторые HTML/HDF5 и checkpoint arrays могли остаться локальными ignored
файлами. Сначала `stat/list` фактических путей. Наличие ссылки в README и зелёного
CI не доказывает наличие локального production dataset.

## Новая организация

Предпочтительно отдельный пакет `src/sloshing/pinned_wetting/` для нового state,
contact handling, detection и output adapter. Исторические directories read-only.
Разрешены необходимые целевые изменения кода, но результаты разных solver source
не должны смешиваться в одной заявленной trajectory.

## Что не делать

Не продолжать STEP3A11/12 I/O investigation. Не пересчитывать все прошлые tests
на production grids. Не навязывать target solver timestep 4.66e-12 из вспомогательного
stiff-bump test. Не говорить «на воде» для nu=0.01. Не обещать molecular resolution.
