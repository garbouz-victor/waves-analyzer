# Вязкий sloshing с закреплёнными контактными точками

Отдельный исследовательский проект на Python. Текущий этап — **STEP 2**:
научная визуализация уже рассчитанного fine dataset, без изменения PDE.
Созданы пять MP4, storyboard и автономный Plotly explorer.
Интерпретация: **CONDITIONALLY QUALIFIED FOR BULK MODEL VISUALIZATION**.
Никакого заранее заданного движения или
экспоненциального затухания в solver нет.

Архив первого этапа с его измерениями и warnings:
[STEP1_REPORT.md](STEP1_REPORT.md).

Независимая проверка и обоснование выбора integrator/dt:
[STEP1_5_REPORT.md](STEP1_5_REPORT.md).

Qualification исходного dataset: [STEP1_6_REPORT.md](STEP1_6_REPORT.md).
Новая область интерпретации и результаты rendering: [STEP2_REPORT.md](STEP2_REPORT.md).

## STEP 2 — scientific bulk visualization

> The animation visualizes the validated bulk behavior of the current
> linear pinned-contact model. The shaded contact zones are deliberately
> excluded from physical interpretation.

STEP 1.6 отклонил **глобальную** small-slope qualification: fine contact slope
достигает 0.3734. Этот исторический verdict не изменён. STEP 2 показывает
согласованные bulk-поля и явно помечает неразрешённые контакты, а не объявляет
весь интерфейс физически разрешённым. Fine — reference, не точное continuum solution.

Источник: `validation_results/animation_qualification/runs/fine-full.h5`:
alpha=0.02°, nu=0.01 m²/s, fine 120×210, SDIRK2, dt=0.00125 s,
snapshot_dt=0.005 s, 0…5 s. Medium используется только для provenance/uncertainty.
Новый smaller-alpha PDE run не нужен для выбранной **bulk** интерпретации.

Скорость вычислена на неподвижной области z≤0; показанная поверхность —
first-order interface displacement, не moving-mesh CFD. Eta показывается
без увеличения: количественная панель в mm, линия в flow panel ×1.
Линейное лагранжево смещение трассеров увеличено **×181**, одинаково по x и z.
Это диаграмма смещения, не физическая экскурсия на сантиметры. Скорость/omega
имеют физические единицы и постоянные цветовые шкалы. Стрелки масштабируются
одним фиксированным коэффициентом с physical velocity key; затухание не скрыто.

Зоны из [animation_scope.json](configs/animation_scope.json): bulk |x|≤0.90,
transition 0.90<|x|≤0.98, unqualified contact |x|>0.98. Это соглашение
визуализации, не резкая граница физической применимости. Серые штрихованные
полосы остаются в explained/clean видео. Ни eta, ни omega не сглаживаются.

### Rendering commands and inspection gates

Нужны `ffmpeg` и `ffprobe` в PATH. Matplotlib, Pillow и Plotly включены в
`requirements.txt`; для package installation: `python -m pip install -e '.[test,render]'`.
Для headless окружения можно задать `MPLCONFIGDIR=/tmp/sloshing-step2-mpl`.
Все команды запускаются из `sloshing_visualization`:

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python scripts/render_step2.py \
  --dataset validation_results/animation_qualification/runs/fine-full.h5 \
  --medium validation_results/animation_qualification/runs/medium-full.h5 \
  --output output/step2 --preview

# Просмотреть 10 PNG в frames_preview/, включая шесть measured key phases.
python scripts/render_step2.py --accept-preview --preview-video
# Просмотреть preview_0_1s.mp4. Только после проверки:
python scripts/render_step2.py --accept-preview-video --all

# Полное декодирование пяти видео, cadence и FEM/source проверки:
python scripts/inspect_step2.py
# Необязательный настоящий browser smoke test: Node >=22 + Chrome.
node scripts/check_step2_explorer.mjs output/step2

python -m pytest -q --junitxml=output/step2/test_runs/unit_final.xml
python -m pytest -q -m validation --junitxml=output/step2/test_runs/validation_final.xml
python scripts/audit_step2_artifacts.py
```

Подрежимы: `--main`, `--clean`, `--vorticity`, `--surface`, `--tracers`,
`--storyboard`, `--explorer`, `--all`. Последний не обходит preview gates.
Если используется другой `--output`, его нужно передавать на каждом шаге.
Принятие preview — явное подтверждение просмотра, не автоматическая оценка красоты.

Renderer **никогда не запускает solver**. При отсутствующих HDF5 остановится:
восстановление выполняется отдельно через STEP 1.6 pipeline, например
`python scripts/qualify_animation_dataset.py --phase all`.
COMPLETE compatible HDF5/display/MP4 кэши переиспользуются. Изменённый fingerprint,
running/failed или несовместимый кэш вызывает ошибку; используйте новый output
либо сначала явно разберите старые артефакты. Полный HDF5 в RAM не загружается.
Карты FEM интерполяции строятся один раз; следующий этап читает display cache по кадру.

Результаты:
[explained](output/step2/bulk_flow_explained.mp4),
[clean](output/step2/bulk_flow_clean.mp4),
[vorticity](output/step2/bulk_vorticity.mp4),
[true-scale surface](output/step2/free_surface_true_scale.mp4),
[tracer comparison](output/step2/tracer_model_comparison.mp4),
[storyboard PDF](output/step2/key_phases.pdf),
[explorer](output/step2/explorer.html),
[metadata](output/step2/animation_metadata.json),
[content fingerprints / audit](output/step2/artifact_audit.json).
MP4, HTML и большой display cache локальны и ignored: после clone эти ссылки
на regenerated artifacts работают после rendering. Metadata, отчёт и несколько
preview PNG/PDF включены в Git. Все MP4: 1920×1080, 50 fps, 1001 snapshot,
20.02 s, замедление физического времени ×4.

Технические API: [Matplotlib FFMpegWriter](https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FFMpegWriter.html),
[Plotly heatmap](https://plotly.com/javascript/reference/heatmap/).
У display heatmap заданы явные границы ячеек, поэтому ни Matplotlib, ни Plotly
не продолжают верхнюю полуячейку в z>0. Explorer содержит 101 реальный snapshot,
работает offline и не является источником scientific error norms.

## Запуск и воспроизводимость

Из каталога `sloshing_visualization`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q

# Короткая проверка solver, не production simulation:
python scripts/run_simulation.py --mesh coarse --t-end 0.1

# Исторический STEP 1 short run (не qualification для animation):
python scripts/run_simulation.py --mesh coarse --t-end 1

# Исторические defaults для воспроизводимости; near-contact slope может
# нарушить линейное приближение. Это НЕ рекомендация для STEP 2 animation:
python scripts/run_simulation.py \
  --a 1 --d 10 --alpha-deg 2 --nu 0.01 \
  --t-end 5 --dt 0.0025 --mesh medium

# Конфигурация + явное переопределение CLI:
python scripts/run_simulation.py --config configs/default.json --nu-preset 0.1
```

Можно установить пакет: `python -m pip install -e '.[test]'`, тогда доступен
`sloshing-run`. `--help` показывает параметры. `--nx` и `--nz` задают число
интервалов независимо от preset. `snapshot_dt` и `t_end` должны быть кратны
`dt`: скрытых изменений временной сетки нет. Последний момент сохраняется,
даже если не совпал с очередным интервалом snapshots.

Все параметры находятся в `SimulationConfig`, включая допуски физических
проверок. `requirements-tested.txt` фиксирует проверенное окружение STEP 1.5/1.6;
`requirements.txt` задаёт совместимые диапазоны для установки, в том числе
на более новых Python. Локальное окружение STEP 1 использует доступные
системные научные библиотеки через `--system-site-packages`; система не
изменялась. Рекомендуемая новая установка выше изолирована полностью.

Имя HDF5 включает вязкость, угол, разрешение и hash всей конфигурации.
Файл хранит JSON параметров, версии зависимостей и время создания. Существующий
run не перезаписывается: для повтора укажите новый `--output results/repeat.h5`.

## Physical problem

Неподвижный прямоугольник `-a ≤ x ≤ a`, `-d ≤ z ≤ 0`, ось z вверх.
По умолчанию `a=1 m`, `d=10 m`, `g=9.81 m/s²`, `alpha=2°`,
`nu=0.01 m²/s`. В момент освобождения:

\[
\eta(x,0)=x\tan\alpha,\qquad \boldsymbol v(x,z,0)=0.
\]

После снятия ускорения объёмной горизонтальной силы нет. Наклонённая
поверхность задаёт неравновесное нормальное напряжение. Исследовательские
вязкости `0.001`, `0.01`, `0.1`, `1.0 m²/s` не называются вязкостью воды.

Уравнения после вычитания гидростатического давления:

\[
\partial_t\boldsymbol v=-\nabla q+\nu\Delta\boldsymbol v,
\qquad \nabla\cdot\boldsymbol v=0,\qquad \boldsymbol v=(u,w).
\]

Адвекция скорости отброшена согласно линейной теории; вязкость непосредственно
входит в PDE при любом положительном nu. Потенциальной аппроксимации нет.

## Boundary conditions and pinned contacts

На обеих боковых стенках и дне **обе** компоненты P2 скорости обращаются
в нуль сильным исключением степеней свободы, включая углы. На `z=0`:

\[
\partial_t\eta=w,\qquad
u_z+w_x=0,\qquad q-2\nu w_z=g\eta,\qquad\sigma=0.
\]

Из no-slip непосредственно следует `w(±a,0)=0`, поэтому концы поверхности
равны начальным значениям при каждом шаге. Они не закрепляются повторным
«исправлением» eta после расчёта: соответствующие строки оператора следа
нулевые после исключения пристеночной скорости.

> Pinned contact line: consequence of strict no-slip model

Условия напряжений естественные, в слабом смысле. Они не навязываются
поточечно в углах: там встречаются разные типы граничных условий и возможны
сингулярности. Никакого contact-angle условия или сглаживания eta не добавлено.

## FEM discretization and pressure

На треугольниках используются непрерывная векторная P2 скорость,
непрерывное P1 давление (Taylor–Hood) и непрерывная P2 поверхность — **тот же
след**, что у вертикальной скорости. Сетка градуирована, квадратура порядка 6.

Слабая форма для тестовой скорости \(\boldsymbol\phi\), нулевой на стенках:

\[
(\partial_t\boldsymbol v,\boldsymbol\phi)
+2\nu(D\boldsymbol v,D\boldsymbol\phi)
-(q,\nabla\cdot\boldsymbol\phi)
+g\langle\eta,\phi_z\rangle_\Gamma=0,
\qquad (r,\nabla\cdot\boldsymbol v)=0.
\]

Здесь \(D\boldsymbol v=(\nabla\boldsymbol v+\nabla\boldsymbol v^T)/2\),
а \(\phi_z\) означает вертикальную компоненту тестовой скорости.
Полный тензор деформаций даёт оба условия свободного напряжения естественно.

**Не нужно дополнительно фиксировать среднее давление.** При заданной eta
его константа определяется нормальным напряжением на открытой поверхности.
Постоянное давление даёт ненулевой функционал
`-c ∫Γ phi_z dx`, поскольку допустимая тестовая скорость там не обязана
иметь нулевой поток. Поэтому у смешанной матрицы нет обычной константной
неопределённости полностью закрытого Stokes. Удаление одного уравнения
давления здесь могло бы нарушить баланс объёма. Тест постоянной поверхности
независимо проверяет `v=0, q=g*eta` и отсутствие ошибочной фиксации давления.

Taylor–Hood обеспечивает несжимаемость относительно пространства P1 тестов,
**не поточечно**. Поэтому и алгебраическая/слабая невязка, и настоящая
`||div v||_L2` сохраняются отдельно; последняя должна изучаться при сгущении.
Она вычисляется из производных FEM по объёмной квадратуре.

Вспомогательная документация: [scikit-fem API: Basis, FacetBasis,
sym_grad](https://scikit-fem.readthedocs.io/en/stable/api.html) и
[операции с производными и проекциями](https://scikit-fem.readthedocs.io/en/stable/howto.html).

## Time integration and discrete energy identity

Implicit midpoint, `dt=0.0025 s`. Обозначим матрицы массы скорости M,
вязкости K, дивергенции B, массы поверхности S; R — вертикальный след,
`C=RᵀS`. После исключения фиксированных скоростей:

\[
M\frac{v^{n+1}-v^n}{\Delta t}+Kv^{n+1/2}
-B^Tq^{n+1/2}+gC\eta^{n+1/2}=0,
\quad Bv^{n+1}=0,
\quad\eta^{n+1}-\eta^n=\Delta t Rv^{n+1/2}.
\]

Подстановка последнего уравнения даёт седловую систему с верхним блоком
`M/dt + K/2 + dt*g*RᵀSR/4`. Sparse LU переиспользуется при всех шагах.
Поверхность связана со скоростью неявно; схемы с явной гравитационной силой нет.

\[
E^n=\tfrac12(v^n)^TMv^n+\tfrac g2(\eta^n)^TS\eta^n,
\qquad
E^{n+1}-E^n=-\Delta t(v^{n+1/2})^TKv^{n+1/2}.
\]

Это проверяется на **каждом внутреннем шаге**, а накопленные и максимальные
невязки сохраняются в snapshots. Энергия и мощность диссипации приводятся
на единицу плотности и единицу поперечного размаха: `m⁴/s²`, `m⁴/s³`.
Для размерной энергии на метр размаха надо умножить на плотность.
Endpoint dissipation в snapshot не подменяет midpoint dissipation в балансе.

Давление в HDF5 относится к **моменту snapshot**, включая t=0. Оно независимо
восстанавливается из мгновенной системы для ускорения с `B dv/dt=0`.
Начальное давление не задаётся нулём. Дополнительная LU факторизация
этой системы также переиспользуется.

Midpoint A-устойчива, но не L-устойчива: очень быстрые вязкие сеточные моды
могут требовать уменьшения dt ради точности. Убывание энергии само по себе
не доказывает точность разрешения пограничного слоя.

## Mesh grading

| Preset | Nx × Nz интервалов | Треугольников |
|---|---:|---:|
| coarse | 40 × 70 | 5 600 |
| medium | 80 × 140 | 22 400 |
| fine | 120 × 210 | 50 400 |

Горизонтальная координата `x=a*tanh(beta*s)/tanh(beta)`, `s∈[-1,1]`,
`beta=1.8`; расстояние вниз от поверхности
`y=d*expm1(gamma*r)/expm1(gamma)`, `r∈[0,1]`, `gamma=4`.
Это сгущает сетку у обеих стенок, поверхности и верхних углов; глубинные
интервалы крупнее. Нулевой grading даёт равномерную координату. Направления
диагоналей зеркальны относительно x=0, чтобы не вносить искусственную
асимметрию между контактами. При нечётном Nx центральная колонка разбивается
на четыре треугольника через центр каждой ячейки (добавляются 2*Nz треугольников).

Tensor mesh также сохраняет узкие пристеночные ячейки на глубине: это
простая первая версия, без адаптивного удаления глубинных вертикальных
линий. Sparse direct solver на fine может требовать существенную память.
В STEP 1.5 выполнен short contact/bulk study на всех трёх presets до 0.1 s.
Это ещё не пространственная сертификация всех полей на интервале 0…5 s.
Малые сетки unit tests не следует считать production-разрешением.

## Diagnostics and failure policy

Каждый сохранённый snapshot имеет строку в
`results/<run>/diagnostics.csv` и те же поля внутри HDF5:

| Проверка | Величина / стандартный допуск |
|---|---|
| No-slip | Шесть точных max абсолютной P2 скорости на стенках, включая экстремумы внутри рёбер, `1e-12 m/s` |
| Слабая несжимаемость | `sqrt((Bv)ᵀ M_p⁻¹ (Bv)) < 1e-9 m/s` |
| Сильная несжимаемость | `||div v||_L2`, и `||div v||/||grad v||`; warning при относительной > 0.15 |
| Объём | `abs(∫eta−∫eta_initial) < 1e-10 m²` |
| Контакты | Максимальное изменение endpoints `< 1e-12 m` |
| Энергия | Относительная ошибка накопленного баланса `< 1e-8`; максимальная ошибка шага и рост также проверяются |
| Наклон | Точный максимум производной каждого P2 участка, левый/правый наклон у стены |

Нулевые все коэффициенты P2 на ребре означают нулевую скорость **на всём
ребре**, а не только в узлах. Unit test дополнительно проверяет сотни
независимых точек стенок. Объём тестируется также независимой точной
квадратурой Симпсона по P2 поверхности. Диссипация независимо сверяется
с интегралом `2 nu D:D` из производных FEM.

Предупреждение об относительной div дополнительно требует абсолютную
`divergence_l2 > weak_divergence_tolerance`: отношение ошибок округления
в практически неподвижном состоянии не считается недоразрешением. Обе
нормы и `velocity_gradient_l2` всё равно сохраняются без обнуления.

Ошибки solver/нефинитные значения/нарушения инвариантов вызывают исключения.
Частичный HDF5 получает `status=failed` и текст ошибки, завершённый —
`status=complete`. Файл с `failed` нельзя выдавать за завершённую simulation.
Warnings не подавляются и не исправляются изменением картинки.

### Small-slope convention

Eta уже размерная: `eta=x*tan(alpha)` на старте. Физически малым должно
быть непосредственно `max |eta_x|`. Если `max_slope` превышает настраиваемый
`small_slope_limit=0.3`, выводится:

> Local small-slope assumption is breaking down near the wall.

Этот порог — диагностический выбор, не резкая физическая граница. Solver
по-прежнему решает заявленные линейные уравнения; физическая интерпретация
локальной области после warning ограничена. В будущем renderer обязан
показывать warning, контакты и реальные наклоны без сглаживания.

## Saved solution and architecture

```text
src/sloshing/
  config.py              параметры и воспроизводимое имя run
  mesh.py                градуированная треугольная сетка
  fem_spaces.py          P2/P1, P2 след, матрицы weak formulation
  time_integrator.py     monolithic midpoint/SDIRK2 и восстановление q(t)
  problems.py            нулевые production loads / явный validation interface
  solver.py              шаги и проверенные snapshots
  diagnostics.py         физические нормы, баланс, warnings
  storage.py             HDF5 и CSV без зависимости от renderer
  cli.py                 config/CLI
  postprocess/
    vorticity.py         dw/dx − du/dz из FEM
    regular_grid.py      интерполяция полей в xv,zv
  validation/            аналитический MMS, stiff/time/contact исследования
scripts/run_simulation.py
tests/                   физические и численные проверки
configs/default.json
results/                 рассчитанные поля
output/step2/            научные видео, metadata, storyboard, explorer
validation_results/      CSV/JSON/PDF/PNG проверок; локальный HDF5 кэш
```

HDF5 содержит `times`; группу `mesh` с координатами, connectivity, dof
coordinates и element dofs для реконструкции P2/P1/DG1; `fields/u,w,q,eta,omega`;
`diagnostics/*`; `probes/coordinates,u,w`; `visualization/xv,zv,u,w,q,omega`.
Формат FEM полей — `(time, dof)`, regular grid — `(time, z, x)`.
`mesh/surface_x` — возрастающие координаты P2 eta, включая середины рёбер
и закреплённые endpoints. Будущий renderer должен оценивать **P2** профиль,
а не заменять его сглаживающим сплайном. Интерполяция сетки отключается
флагом `--no-visualization-grid` для облегчённых тестов.

Завихренность `omega=dw/dx−du/dz` — точная поэлементная линейная производная
квадратичной скорости. Она хранится в разрывном P1 (DG1), через точную L2
проекцию; усреднения через границы элементов нет. На границе элементов
regular grid выбирает одну из односторонних величин; для диагностики надо
использовать FEM данные, не различать пиксели. Регулярная сетка 81×241 удобна
для будущего renderer, но сама по себе недостаточна для самых тонких слоёв:
contact zoom должен опрашивать FEM с дополнительным разрешением.

## How to watch the animation

Для текущей STEP 2 анимации:

1. Сначала смотрите на TRUE-SCALE eta в mm: справа выше, затем bulk профиль
   проходит через почти горизонтальное положение; endpoints не двигаются.
2. Сопоставьте velocity arrows: после release вверху поток влево, слева вверх,
   справа вниз. Цвет speed и длины стрелок имеют постоянные шкалы.
3. Сравните с K: первый максимум t=0.395 s близок к modal crossing t=0.4161 s,
   но это не один и тот же момент и не заданное равенство.
4. Следите за слоями linear tracers ×181 и full-depth профилем Q(z):
   движение резко ослабевает с глубиной.
5. После overshoot следите за reversal стрелок и последующим затуханием.
6. В отдельном vorticity video смотрите на изменение знака и пристеночные слои.
7. Серые contact strips не интерпретируйте как разрешённую contact-line physics.

При небольшой вязкости возможны overshoot и смена направления циркуляции;
при большой — сильное подавление колебаний. В глубоком сосуде движение
может быть сосредоточено сверху. Наличие и глубина обратного течения,
время максимумов энергии и картина циркуляции требуют проверки полей:
несжимаемость сама по себе не предписывает знак горизонтальной скорости
на каждой глубине. Они не закладываются в визуализацию как сценарий.

Default STEP 2 tracer: `X0 + integral v(X0,t) dt`, скорость оценивается
в фиксированной начальной точке. Это first-order Lagrangian displacement.
Отдельный comparison использует RK4 pathlines в линейном Eulerian field;
их отличие higher-order и не доказывает nonlinear drift/mixing. Цвет начальной
глубины постоянен; invalid trajectories не отражаются и не clip'аются.

> Streamlines = instantaneous velocity field
>
> Pathlines = trajectories integrated in the supplied unsteady velocity field;
> here this does not validate nonlinear transport in the physical experiment.

## What this animation can and cannot tell us

Эта модель позволяет исследовать bulk sloshing, поле скорости, вязкое
затухание, генерацию завихренности, пристеночные слои и последствия pinned
contact lines. Это двумерная линейная модель на неподвижной области.

Она **не решает** deposition of a film on a vertical wall, dewetting,
moving contact line или final wall drainage. Строгий no-slip закрепляет
контакты; single-valued graph `z=eta(x,t)` не описывает вертикальную плёнку.
Для движения контактной линии потребуется отдельная wetting/contact-line
model. Фиктивная плёнка в результаты не добавляется.

При sigma=0 и неравных закреплённых высотах гладкое конечное состояние покоя
с горизонтальной поверхностью не удовлетворяет обоим endpoints. Это делает
предел больших времён и локальные градиенты у контактов особенно чувствительными
к сетке. Нельзя приписывать остаточную пристеночную структуру физической плёнке
или считать проверку убывания энергии достаточной проверкой её сходимости.

## STEP 1.5: independent validation commands

```bash
python -m pytest -q
python -m pytest -q -m validation
python scripts/generate_manufactured_solution.py
python scripts/manufactured_convergence.py
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/time_stiffness_study.py --workers 2
python scripts/contact_mesh_study.py
```

Обычный pytest включает symbolic MMS и короткий FEM MMS; marker `validation`
включает дорогую пространственную сходимость и проверку загрязнения dt.
`time_stiffness_study.py --unit-only` ограничивает исследование скалярной
задачей и tiny FEM-модой. Полная команда выполняет 18 medium run до 0.5 s.
CSV, JSON, научные PDF/PNG находятся в `validation_results/`. Кэш HDF5 содержит
FEM коэффициенты для повторного анализа и не включается в git. Завершённые
совместимые runs переиспользуются; failed/incompatible cache вызывает ошибку.
Для независимого повторения укажите новый `--output`.
`--workers 2` запускает независимые вязкости в отдельных процессах; это
не меняет метод интегрирования. `--physical-only` переиспользует уже
сделанные scalar/FEM-mode studies и запускает только physical refinement.
Необрабатываемое завершение процесса может оставить `running` или повреждённый
HDF5; такие файлы не считаются завершёнными и автоматически не перезаписываются.

### Analytic manufactured solution

MMS выполняется на `a=d=1`, равномерной сетке; production глубина остаётся 10 m.
`psi=exp(-t)*(1-x²)²*(z+1)²`, `u=psi_z`, `w=-psi_x`,
`eta=exp(-t)*d((1-x²)²)/dx`, `q=exp(-t)*x*(z+1)`.
SymPy независимо выводит `f=v_t+grad(q)-nu*Laplacian(v)` и
`r=Tn+g*eta*n`. Они входят как `(f,phi)+<r,phi>` в правую часть той же weak
formulation, включая обе компоненты traction. Exact fields не зависят от FEM
матриц. [SymPy lambdify](https://docs.sympy.org/latest/modules/utilities/lambdify.html)
преобразует символические выражения в NumPy функции.

ProductionProblem возвращает нулевые силы; CLI не имеет переключателя на MMS.
Validation явно создаёт ManufacturedProblem. Initial velocity — constrained
L2 projection аналитической скорости, initial eta — P2 L2 projection с точными
endpoints. Это исключает несовместимую начальную div в DAE. Для нагрузки
используется квадратура порядка 10, для ошибок — 12. Ошибки считаются
относительно аналитических функций без сдвига константы pressure и без
visualization grid; velocity H1 в таблицах означает градиентную полунорму.

### Monolithic SDIRK2 and its energy budget

`--integrator midpoint|sdirk2`, default `midpoint` сохранён. Обозначим
`gamma=1-1/sqrt(2)`, `s=gamma*dt`, `b=(1-gamma,gamma)`. Для стадии i:

```
v_base   = v_n   + dt * sum_{j<i} a_ij k_v,j
eta_base = eta_n + dt * sum_{j<i} a_ij R V_j

(M/s + K + s*g*C*R) V_i - B^T Q_i = M*v_base/s - g*C*eta_base + F(t_n+c_i*dt)
B V_i = 0
eta_i = eta_base + s*R V_i
k_v,i = (V_i-v_base)/s
```

`c=(gamma,1)`, `a_21=1-gamma`, конечное состояние равно второй стадии.
Используется одна LU факторизация на обе стадии и все шаги. Никаких
проекций div после шага или явного обновления eta нет. Midpoint использует
силу в середине шага; SDIRK2 — в соответствующие моменты стадий.

Для общего RK точный дискретный баланс квадратичной энергии:

```
E_{n+1} - E_n + dt*sum_i b_i V_i^T K V_i - dt*sum_i b_i V_i^T F_i + Q_RK = 0
Q_RK = dt²/2 * sum_ij (b_i*a_ij+b_j*a_ji-b_i*b_j) <k_i,k_j>_G
G = diag(M, g*S)
```

У midpoint `Q_RK=0`. У данного SDIRK2 матрица в скобках равна
`diag(-gamma²,+gamma²)`. Поэтому поправка **знаковая**, а метод не является
algebraically stable. Это не отменяет его A-/L-stability для линейной задачи.
Поправка учитывается отдельно, не выдаётся за положительное вязкое тепло.
Сохраняются `cumulative_viscous_dissipation`, `cumulative_work`,
`rk_energy_correction`, итоговая невязка баланса; для production дополнительно
проверяется отсутствие роста полной энергии на каждом внутреннем шаге.

Завихренность сравнивается во всём объёме, слое `z>=-1` и полосах `x<=-0.9`,
`x>=0.9`. Интегралы DG1² по обрезанным треугольникам точные. Для выявления
temporal alternation сохраняются нормы и probes на **каждом dt**, включая
точки `(±0.995,-0.005)` около контактов. Показатель `nu*dt/h_min²` основан
на минимальной длине ребра и является только rough dimensional indicator,
не спектральным радиусом FEM.

Единицы новых норм: `||v||_L2` — m²/s, `||grad v||_L2` и `||omega||_L2`
в двумерной области — m/s, `||q||_L2` — m³/s², поверхностная `||eta||_L2`
— m^(3/2), `omega_max` и omega probes — 1/s. Относительные ошибки — доли;
в отчёте проценты явно помечены `%`. MMS — математический benchmark линейного
оператора, его большая аналитическая амплитуда не интерпретируется как
физически допустимый свободный интерфейс.

`configs/validation_alpha_0_2.json` (ранее `linear_safe.json`) задаёт alpha=0.2°, исторический **кандидат** для линейной
физической анимации; это имя не гарантирует малого наклона у контакта во все
моменты времени. `configs/contact_breakdown_demo.json` сохраняет alpha=2°:
полезный формальный пример быстрого нарушения малого наклона у pinned contacts,
а не «неправильный solver run». Измеренные ограничения и рекомендуемые
настройки приведены в STEP1_5_REPORT. При nu=0.01, medium, alpha=0.2°
измерен max slope=0.879 к 0.5 s: даже этот кандидат нарушает локальное
условие малого наклона. Название конфига не является физической гарантией.

### STEP 1.5 measured recommendation (historical, not an animation release)

Для следующего физического расчёта рекомендуются **SDIRK2, dt=0.00125 s**,
medium, snapshots через 0.005 s. У SDIRK2 значительно меньше stiff residue
в omega, чем у midpoint; dt=0.0025 уже хорошо согласуется по нормам полей
на snapshots, но ближайшие к контакту probes в первые миллисекунды требуют
осторожности. Default integrator и default alpha не изменены.

Для физической интерпретации контактов разумный следующий угол — **0.02°**:
по линейному масштабированию измеренного решения это даёт max slope≈0.0879
на medium за 0…0.5 s при nu=0.01. Это вывод из линейности, не отдельный run
и не гарантия для fine или 5 s. Alpha=0.2° остаётся параметром выполненного
validation study и формального исследования локального breakdown.

Подтверждены P2/P1 MMS orders, временная сходимость на medium и уменьшение
bulk/contact L2 differences на трёх сетках при t=0.1. Pointwise corner slope
и omega_max **не** объявляются пространственно сошедшимися. Интервал 0.5…5 s
и nu=0.001 пока не прошли это temporal validation. Подробные числа и полный
перечень тестов находятся в [отчёте](STEP1_5_REPORT.md).

## STEP 1.6 historical qualification

STEP 1.6 использует `configs/animation_candidate.json` и
`configs/animation_candidate_fine.json`: alpha=0.02°, nu=0.01, SDIRK2,
dt=0.00125, snapshot_dt=0.005, t_end=5. Эти configs — кандидаты, не обещание
малого slope на всём интервале. Измерения и решение:
[STEP1_6_REPORT.md](STEP1_6_REPORT.md). На этапе STEP 1.6 production animations отсутствовали.

Результат реально выполненных medium/fine 0…5 s: **NOT QUALIFIED FOR
PHYSICAL ANIMATION** для alpha=0.02°. На fine max slope=0.373408,
первый сохранённый snapshot выше 0.3 — t=2.365 s. Это локальное нарушение
linear-small-slope policy у pinned contact, не instability solver. Bulk eta
согласуется в пределах 0.0760%, velocity — 0.0793%, wall omega L2 — 5.731%
от пиковых fine norms. Максимальная bulk path separation — 0.403 µm.
No-slip, volume, energy и symmetry проходят; strong divergence уменьшается
на fine, но остаётся отдельной оговоркой для medium.

В STEP 1.6 предложен **не рассчитанный** conservative candidate: alpha=0.0025°, fine,
nu=0.01, SDIRK2, dt=0.00125, snapshot_dt=0.005, t_end=5. По tan scaling
измеренного fine решения ожидался max slope≈0.04668. STEP 2 не выполняет этот
run: уменьшение amplitude не разрешает continuum corner и не улучшает
relative field errors. Вместо глобальной qualification принята явно ограниченная
bulk interpretation, описанная выше. Исторический отчёт/summary сохранены без
изменений; текущие configs остаются выполненными 0.02° cases.

### STEP 1.6 commands and gates

```bash
python scripts/qualify_animation_dataset.py --phase medium-short
python scripts/qualify_animation_dataset.py --phase medium-full
python scripts/qualify_animation_dataset.py --phase fine-short
python scripts/qualify_animation_dataset.py --phase fine-full
python scripts/qualify_animation_dataset.py --phase compare
python scripts/qualify_animation_dataset.py --phase particles
python scripts/qualify_animation_dataset.py --phase summary
# All gates sequentially, or resume using COMPLETE compatible caches:
python scripts/qualify_animation_dataset.py --phase all
```

Успешный exit команды означает завершение исследования; физический verdict
нужно читать в `summary.json` (`qualified` / `conditional` / `failed`).

Для каждой команды рекомендуется окружение `OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`. Default output:
`validation_results/animation_qualification`; другой `--output` даёт новый
независимый расчёт. Full run продолжает COMPLETE short prefix в **новом**
файле, сохраняя исходный prefix. State restart включает cumulative viscous
loss, RK correction и энергетические максимумы. Running/failed/truncated
HDF5 не переиспользуется. Это фазовый restart, не автоматическое возобновление
произвольно оборванной записи. Источник comparison проверяется по config,
размеру и времени изменения completed dataset. При первом чтении каждой
ревизии файла проверяются все FEM snapshots и probes на конечность и полноту;
не только конечный state. Новые chunked datasets используют Fletcher32.
Старые completed prefixes не переписываются ради добавления checksum.

HDF5 содержит P2 u/w, P1 q, P2 eta, DG1 omega, mesh/dof maps, probes,
diagnostics и restart accounting. Visualization grid не сохраняется.
`runs/`, `raw/` и NPZ trajectories ignored; в Git — компактные summaries
и научные PDF/PNG. Units: regional slope — безразмерный; volume L2 velocity —
m²/s; omega/div L2 — m/s; surface L2 — m^(3/2); horizontal velocity L2(dx) —
m^(3/2)/s; particle separations — m.

Cross-mesh integrals используют точное геометрическое разбиение по рёбрам
обеих сеток, затем квадратуру порядка 4 на общих треугольниках. Это точно
интегрирует квадрат разности P2 velocity и DG1 omega. Regional P2 slopes
вычисляются по концам соответствующих отрезков, включая границы областей
внутри FEM edge. Угловые производные не сглаживаются.

Контрольные частицы: 28 bulk + 12 near-wall seeds. Spatial interpolation —
реальные P2 elements; time interpolation — линейная между snapshots; RK4
останавливает шаг на каждом временном узле. Проверяются RK step refinement,
snapshot_dt=0.005 vs 0.0025 при одном PDE dt и medium vs fine. Выход над
P2 surface/reference domain делает путь invalid навсегда, без отражений;
first-invalid time определяется с точностью RK stage, не event root solver.

Qualification policy отделена от solver tolerances. Slope<=0.1 — conservative,
0.1…0.3 — conditional, >0.3 — failed linear-small-slope qualification.
Инженерные ориентиры для интегрального сравнения: bulk eta 1%, velocity 2%,
wall omega 10%; bulk path mesh difference 0.005 m, snapshot difference 10 µm,
ODE-step difference 1 µm.
Они заданы до получения medium/fine comparison, не являются теоремами
об ошибке и публикуются вместе с фактическими значениями. Significant
strong-divergence caution: >5% более чем на 10% активного времени; активные
snapshots имеют gradient norm >1% пикового. Уменьшение alpha снижает
абсолютный slope, но не улучшает относительные пространственные ошибки
или relative strong divergence.

GitHub Actions выполняет только `python -m pytest -q` на Python 3.11,
не длинные qualification runs. Конфигурация основана на официальных
[setup-python](https://github.com/actions/setup-python) и
[checkout](https://github.com/actions/checkout); локальные научные измерения
пока сделаны в указанном в requirements-tested.txt окружении Python 3.9.

STEP 2 является отдельной итерацией после явного решения qualification.
Прохождение unit tests само по себе не разрешает физические утверждения
об анимации. Measured medium/fine fields и проверка траекторий выполнены,
но глобальный small-slope gate не пройден. Неразрешённый contact corner
нельзя выдавать за точную физику стенки.
Production renderer, MP4 и Plotly не входят в STEP 1.6.
