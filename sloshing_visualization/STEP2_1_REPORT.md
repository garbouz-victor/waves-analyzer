# STEP 2.1 — clarify no-slip boundary visualization

## 1. User-visible problem

До изменения renderer проведён аудит commit
`4172787060cdf38ed590eaaca393dc6fb6e8c0e7`. На левой стенке в
bulk_flow_explained кажется возможным «ползучее» движение. Это неоднозначность
визуального расположения стрелок, не установленное нарушение FEM no-slip.

## 2. Why the old animation looked like slip

`data.py` строит quiver x=linspace(-1,1,25). В x=-1 скорость нулевая,
а `panels.py` подавляет minimum-length dot (`minlength=0`): нулевой вектор
невидим. Ближайшая ненулевая колонка — x=-0.9166667, **8.333 cm внутри жидкости**.
Фиксированный коэффициент 102.334568 s превращает 0.5 mm/s в 0.0511673 m
длины рисунка. Tail/seed не выделен; glyph может восприниматься как скорость
стены. Полная серая вертикальная полоса contact warning дополнительно смешивает
неразрешённый угол с корректно заданной solid boundary.

Incorrect visual inference: “Arrow near wall appears to indicate wall motion.”

Actual numerical interpretation: “The arrow is seeded inside the fluid.
The exact wall FEM velocity is zero.”

**Contact-line singularity != violation of no-slip.**

Первые два раздела записаны **до изменения renderer**. Последующие разделы
содержат измерения actual source/cache, а не гипотезу о виновности solver.

## 3. Verification that FEM no-slip is actually satisfied

Никакая production FEM simulation не запускалась. Исходные PDE, boundary
conditions, коэффициенты FEM и физические параметры не редактировались.
Использован существующий `fine-full.h5` из STEP 1.6: α=0.02°, ν=0.01 m²/s,
a=1 m, d=10 m, fine 120×210, SDIRK2, dt=0.00125 s,
snapshot_dt=0.005 s, 0…5 s, 1001 snapshot. Medium — только provenance.

Полный SHA-256 fine HDF5 после postprocessing:

```text
fa299d74c9487a40a5fac55966c1288c950bc0a46116f80e7da0063f794cf09a
```

Это исходный hash, **без изменения единого байта**. Численный dataset создан
STEP 1.6 commit `dc51e33f…`; базовый renderer commit этой итерации —
`4172787060cdf38ed590eaaca393dc6fb6e8c0e7`. `source_commit_sha` в новом
animation metadata означает этот postprocessing base, не новый расчёт PDE.

| Проверка всех 1001 snapshots | Результат |
|---|---:|
| max wall-left speed, m/s | 0 |
| max wall-right speed, m/s | 0 |
| Независимая проверка cache profiles относительно source P2 | точное совпадение |
| Максимальное расхождение отражённых speed profiles, m/s | 3.1767123654×10⁻¹⁷ |
| Относительная profile asymmetry | 2.7531181×10⁻¹⁴ |

В `display_data.h5` добавлены массивы `wall_left_speed_max(t)` и
`wall_right_speed_max(t)`. Это не значения сглаженной картинки. P2 evaluation
опрашивает все 421 trace nodes каждой стены. На каждом P2 сегменте максимум
нормы векторного полинома находится по концам и стационарным точкам quartic
|v|² (кубическое уравнение для производной). Для этой constrained trace все
nodal coefficients равны нулю, следовательно весь P2 trace точно нулевой.
Функция также протестирована на ненулевом полиноме с максимумом между узлами.

`scripts/inspect_no_slip.py` независимо проверил source walls для всех кадров,
профили для восьми времён, symmetry и geometry за **74.82 s**.
Повторная итоговая проверка с привязкой к renderer signature заняла **261.01 s**
при одновременном кодировании нескольких видео; все числа воспроизвелись.
Обычные physics tests, включая обе стены и дно, сохранены. Unit tests могут
решать tiny задачи; это не пересчёт production dataset.

## 4. Arrow seed geometry before

Старый grid: 25 колонок x∈[−1,1], шаг 0.0833333 m.
В стеночных колонках нет видимых стрелок: velocity=0, `minlength=0`.
Первая **ненулевая** левая колонка: x=−0.9166666667, не x=−1.

До исправления из **реального encoded MP4**, snapshot 79, physical t=0.395 s,
извлечён кадр. Известные data coordinates перенесены на изображение:

| Выбранная стрелка | x, m | z, m |
|---|---:|---:|
| tail / FEM seed | −0.916666667 | −0.150000000 |
| tip | −0.928567457 | −0.067651348 |

Стенка, seed, tail и tip отмечены в
[before_arrow_geometry.png](output/step2/before_arrow_geometry.png),
числа — в [JSON](output/step2/before_arrow_geometry.json).
Не утверждается, что эта стрелка математически пересекала стену: неверным
было визуальное отождествление близкого glyph с boundary velocity.

## 5. Arrow glyph length before

Фиксированный `fixed_arrow_seconds=102.33456782346151 s` — графический
перевод скорости в длину стрелки, **не время физического переноса частицы**.
Physical key 0.5 mm/s соответствует длине 5.11673 cm на диаграмме.

По всем старым cached arrow samples за 0…5 s:

| Величина | Значение |
|---|---:|
| max rendered horizontal excursion, S·|u| | 0.0839121193 m |
| max rendered vector length, S·|v| | 0.0997671112 m |

Эти maxima могут принадлежать разным seeds/временам. Нельзя вычитать
глобальную максимальную длину из положения ближайшей стрелки и объявлять
пересечение: именно поэтому проверена фактическая glyph geometry.

## 6. New arrow geometry

Выбран разрешённый фиксированный interior range **[−0.85,0.85] m**,
25 колонок; 18 z-уровней от −1.45 до −0.15 m.
`pivot="tail"` явно установлен. Ни длины, ни скорости не обрезаются при
приближении к стене; normalization отсутствует; S **не изменён**.
Seeds и scales постоянны во времени.

До рендера вычислена консервативная граница целого glyph:

```text
max |seed_x| + S max |u| + head_half_width
0.85 + 0.0839121193 + 0.00816 = 0.9420721193 m
```

При внешней границе contact strip |x|=0.98 эта оценка оставляет
0.0379278807 m, что больше заданного extra margin 0.02 m.
До solid wall гарантированная консервативная оценка — **0.0579278807 m**.
Учитывается ширина головы, не только точка tip. Никакой adaptive
framewise перестановки seeds нет.

Дополнительно проверены **все вершины всех Matplotlib glyphs во всех 1001
кадрах**, через public `get_paths`, transforms и offsets. Это геометрическая
проверка artist coordinates, не эвристика распознавания пикселей.

| Новая геометрия | Измерение |
|---|---:|
| max S·|u| | 0.0839121193 m |
| max S·|v| на новых seeds | 0.0995201617 m |
| tip x range | [−0.8762711421, 0.8691450380] m |
| весь glyph x range | [−0.8777066955, 0.8711361834] m |
| фактический минимальный glyph–wall зазор | **0.1222933045 m** |
| самый высокий glyph z | −0.0538944036 m |

Неравенство между консервативным и фактическим зазором нормально:
максимальная горизонтальная скорость возникает не у крайнего seed.
Эти data-coordinate расстояния измерены до линии x=±1. Дополнительно учтена
конечная ширина рисунка стены в actual 1920×1080 layout: внешняя белая
обводка занимает до 0.006751543 m внутрь панели. Даже до внешнего **ink edge**
консервативный зазор **0.051176338 m**, фактический минимум по всем кадрам
**0.115541761 m**. Это всё ещё больше extra margin 0.02 m.
Одинаковая гарантия действует в main/clean/vorticity panels. Отдельный
wall zoom использует собственный **фиксированный** S=20 s с подписанным key;
он показывает ближе расположенные seeds без имитации движения стены.

Сравнение actual encoded кадров: [до](output/step2/before_arrow_geometry.png)
и [после](output/step2/after_arrow_geometry.png).
В новом кадре t=0.395 s отмечены tail (−0.85,−0.15) m и
tip (−0.875376514,−0.054276341) m. Исходная картинка «до» сохранена до
изменений; повторное аннотирование archived encoded video со старым layout
воспроизводит её точно, без redraw старого velocity field.

## 7. New wall rendering

Поверх heatmap, quiver, tracers и corner patch рисуются чёрные solid-wall
линии linewidth=2.3 pt, zorder=30 с белым контрастным outline. Это boundary
glyph нулевой толщины в математической модели, не искусственный fluid strip.
На full-depth view также обозначено дно z=−10. В узком left-wall zoom
рисуется только видимая левая стенка: off-screen правая не проецируется
поверх соседних profile panels.

Leader `no-slip: u = w = 0` связан непосредственно со стенкой.
В explained/vorticity версиях показано `FEM wall velocity (actual P2):
left 0, right 0 m/s — no-slip`. Clean сохраняет wall label и геометрию,
но без дополнительного dashboard. Стеночная нулевая скорость не изображается
ненулевой минимальной quiver-точкой.

## 8. Wall-normal velocity profiles

Новая `wall_normal_profiles()` использует actual P2 element mapping,
а не interpolation от display-grid. Слева x=−1+ε, справа x=1−ε.
Глубины −0.05, −0.2, −1 m. ε включает 0, равномерную сетку 0…0.15 m
и дополнительные точки 10⁻⁵,10⁻⁴,5×10⁻⁴,10⁻³,… m.
В cache сохранены отдельно u, w и speed для обеих стенок. Fixed limits
вычислены один раз из полного temporal dataset, без нормирования профилей:

| Profile axis | Постоянные пределы, mm/s |
|---|---:|
| speed | [0, 1.223091402] |
| u | [−0.394445479, 0.394445479] |
| w | [−1.166503453, 1.166503453] |

Примеры actual **speed в mm/s**; во всех строках при ε=0 три значения
ровно 0. Это значения поля, не fitted boundary-layer law.

| t, s | ε, m | z=−0.05 m | z=−0.2 m | z=−1 m |
|---|---|---:|---:|---:|
| 0.100 | 0.00001 | 0.000229421 | 0.000133737 | 0.0000299734 |
| 0.100 | 0.001 | 0.0226617 | 0.0132031 | 0.00295875 |
| 0.100 | 0.01 | 0.201795 | 0.117047 | 0.0262039 |
| 0.395 | 0.00001 | 0.000261652 | 0.000113460 | 0.0000265585 |
| 0.395 | 0.001 | 0.0261134 | 0.0113452 | 0.00265500 |
| 0.395 | 0.01 | 0.254698 | 0.112993 | 0.0263840 |
| 0.810 | 0.00001 | 0.0000248120 | 0.0000944787 | 0.0000300930 |
| 0.810 | 0.001 | 0.00232758 | 0.00930918 | 0.00297133 |
| 0.810 | 0.01 | 0.0168550 | 0.0809001 | 0.0263677 |
| 1.615 | 0.00001 | 0.00000507499 | 0.0000739645 | 0.0000243581 |
| 1.615 | 0.001 | 0.000662368 | 0.00726239 | 0.00240176 |
| 1.615 | 0.01 | 0.0215938 | 0.0609246 | 0.0210219 |
| 4.040 | 0.00001 | 0.00000689435 | 0.0000565600 | 0.0000176241 |
| 4.040 | 0.001 | 0.000589589 | 0.00556124 | 0.00173907 |
| 4.040 | 0.01 | 0.00915242 | 0.0473233 | 0.0153356 |

Видно стремление скорости к нулю при ε→0, зависимость от глубины и фазы.
Не предполагается монотонность speed по всему ε: при reversal profile может
иметь внутренние extrema. Никакая форма v∼ε или power law не подгонялась.
Участок z=−0.05, ε<0.02 лежит в **помеченном corner patch**: это корректные
вычисленные значения текущей дискретной модели, не новая квалификация
continuum contact-corner physics. Более глубокие профили отделяют эту
оговорку от общего математически точного no-slip.

[no_slip_profiles.pdf](output/step2/no_slip_profiles.pdf) и
[PNG](output/step2/no_slip_profiles.png) содержат t=0.100,0.395,0.810,
1.220,1.615,4.040 s — все это реальные snapshots. Малый
[JSON](output/step2/no_slip_profile_samples.json) хранит выборку обоих profiles.

## 9. Left/right symmetry

Линейное antisymmetric release даёт отражательную структуру: u слева/справа
совпадает, w меняет знак, поэтому speed совпадает. Никакая симметризация
cache не выполнялась. Максимальная абсолютная разность speed profiles
слева/справа по всем временам/точкам — **3.1767×10⁻¹⁷ m/s**.
Это roundoff-level discrepancy, не видимое нарушение симметрии.
Synthetic P2 tests сравнивают независимые analytical polynomial values
и `scikit-fem probes`; production regression проверяет компоненты со
знаками отражения, а не только speed.

## 10. Contact corner vs solid-wall no-slip

**Contact-line singularity != violation of no-slip.**

Pinned contact issue: eta endpoint не может двигаться при strict no-slip
и кинематическом условии eta_t=w. Это порождает локальную corner/free-surface
проблему, установленную STEP 1.6. Solid-wall condition u=w=0 при этом
реализована корректно во всём FEM trace.

На SurfacePanel x-strips |x|>0.98 и слабый transition 0.90<|x|≤0.98
остаются уместны: это график поверхности. На 2D Flow/Vorticity/Tracer panels
contact warning заменён угловыми patches:

```text
left:  −1 ≤ x ≤ −0.98, −0.12 ≤ z ≤ 0.04 m
right:  0.98 ≤ x ≤ 1,  −0.12 ≤ z ≤ 0.04 m
```

Высота patch — **visualization convention, не physical singularity thickness**.
Слабое transition shading из основного 2D flow убрано, чтобы не закрывать
пристеночный градиент. Подпись: `Unresolved contact corner; solid walls below
it obey no-slip`. Ни вся стенка, ни весь пристеночный слой не объявляются
нефизичными. Explorer сохраняет предупреждение даже при выключенной заливке.

## 11. Updated videos

Все artifacts находятся в `output/step2/`. Старые generated STEP 2 outputs
сохранены recoverably в ignored `archive_step2_4172787`, не удалены;
production HDF5 никогда не перемещался. Display schema изменилась на
`step2-display-v2-no-slip`; renderer signature включает новые profile/panel
модули. Старые preview acceptance и MP4 больше не совместимы.
Повторная preparation cache заняла **141.22 s**, без запуска solver.

Статический preview повторён и явно просмотрен: шесть modal phases,
clean/vorticity/tracer layouts и специальные enlarged crops
[t=0.100](output/step2/frames_preview/no_slip_t0.1.png),
[0.395](output/step2/frames_preview/no_slip_t0.395.png),
[0.810](output/step2/frames_preview/no_slip_t0.81.png),
[1.615](output/step2/frames_preview/no_slip_t1.615.png).
После статического принятия выполнен новый 0…1 s MP4 (201 frames),
полностью декодирован и просмотрены extracted frames, затем принят video gate.
Один параллельный worker попытался стартовать до записи accepted gate и
правильно остановился **до создания MP4**. Он перезапущен после acceptance;
это защитная проверка последовательности, не numerical failure.

Два длинных render processes получили **exit 143 / SIGTERM**: первый проход
tracer comparison после сообщения о кадре 901, первый no-slip проход после
кадра 329. Причина внешнего завершения не установлена; Python/FEM exception
не сообщался, production solver вообще не работал. Эти MP4 **не приняты**:
их running manifests помечены failed, частичные файлы recoverably перенесены
в ignored `interrupted_step2_1`, после чего рендер повторён из неизменного
validated display cache. Это не замаскированный success и не повтор PDE.
Малый [render_interruptions.json](output/step2/render_interruptions.json)
сохраняет инцидент отдельно от успешных final manifests.

Все **шесть полных MP4 завершены и полностью декодированы** без ошибок.
Каждый: 1920×1080, 50 fps, 1001 frame, duration 20.020 s, физическое время
0…5 s, замедление ×4. Проверены PTS-интервалы 0.02 s; из каждого извлечены
кадры 83 и 1000. Готовое видео не принималось только по наличию файла.

| Обновлённый artifact | Успешный render, s | MiB |
|---|---:|---:|
| [bulk_flow_explained.mp4](output/step2/bulk_flow_explained.mp4) | 1391.27 | 5.351 |
| [bulk_flow_clean.mp4](output/step2/bulk_flow_clean.mp4) | 1305.09 | 4.855 |
| [bulk_vorticity.mp4](output/step2/bulk_vorticity.mp4) | 1233.22 | 4.074 |
| [free_surface_true_scale.mp4](output/step2/free_surface_true_scale.mp4) | 395.47 | 1.149 |
| [tracer_model_comparison.mp4](output/step2/tracer_model_comparison.mp4) | 250.10 | 0.879 |
| [no_slip_boundary_layer.mp4](output/step2/no_slip_boundary_layer.mp4) | 1083.42 | 2.194 |

Времена — фактические wall-clock времена конкретных успешных попыток;
часть процессов выполнялась одновременно, поэтому это не сравнительный
performance benchmark. Суммарный размер шести MP4 **18.501 MiB**. Они,
display HDF5, большой HTML и temporary decoded frames остаются ignored.
Коммитятся source/config/tests, metadata/manifests, отчёт, несколько PNG,
одностраничный profiles PDF и обновлённый девятистраничный key_phases PDF.

Новая подпись renderer/gates:

```text
e4d52c9374ec5665ce559228ea901d09fa9a679f2557fe22011cf814ecc1c1c3
```

[inspection_summary.json](output/step2/inspection_summary.json) содержит
decode/cadence checks и video SHA-256; [no_slip_inspection.json](output/step2/no_slip_inspection.json)
— все-frame FEM/glyph проверки; [artifact_audit.json](output/step2/artifact_audit.json)
— content hashes, неизменность истории, source-code соответствие green CI,
проверку полного набора тестов и запас до внешней обводки стенки.
Команды воспроизведения, повторного принятия gates и проверки — в
[README](README.md#rendering-commands-and-inspection-gates). Renderer solver
не запускает. Для длительных сессий можно запускать отдельные `--main`,
`--clean`, `--vorticity`, `--surface`, `--tracers`, `--no-slip` после принятия
обоих gates; несовместимые или incomplete outputs молча не переиспользуются.

### Tests and CI

`tests/test_contact_mask.py` больше не предполагает, что `get_xy()` возвращает
массив вершин. Маски созданы как Rectangle и проверяются public
`get_x/get_width/get_y/get_height`. Исправление проверено локально на
Matplotlib 3.5.1 и реальным GitHub Actions на **Matplotlib 3.11.1**.

Новый `test_no_slip_visualization.py`: **11 passed**, включая обе стороны
P2 profiles, exact wall evaluation, внутрисегментный максимум P2 trace,
symmetry, seed/tip safety, tail pivot, label, zorder, corner-only geometry,
дно full-depth и fixed speed/u/w profile axes. Полный обычный suite:
**128 passed, 2 deselected**, 54.53 s в финальном локальном запуске.
Отдельный `pytest -q -m validation`: **2 passed, 128 deselected**, 86.45 s.
Итого **130 passed** в двух непересекающихся наборах; 11 из них — новые
no-slip visualization tests. Полный список каждого test id сохраняется в
`output/step2/test_summary.json` из финальных JUnit XML.
Предупреждения прежних coarse physics tests не скрыты: малый slope и strong
divergence на намеренно грубых test meshes; это не warnings production rerun.

Реальный [GitHub Actions run 34057551865](https://github.com/garbouz-victor/waves-analyzer/actions/runs/34057551865)
**green**, 128 passed, 2 deselected, 14 warnings, 6.90 s;
Python 3.11.16 / Matplotlib 3.11.1 / NumPy 2.4.6. Чтобы сначала получить green
CI, а затем создать финальный commit, использован отдельный технический
validation snapshot `5360e12b74bd6b0f27884ab672485d8e4cb30ba1` в ветке
`validation/step2-1-no-slip`, без изменения remote master.
Это не заявление, что старый failed run задним числом стал green.
Полные данные — [ci_summary.json](output/step2/ci_summary.json).

Explorer прошёл настоящий offline Chrome smoke test: slider, play/pause,
fixed speed/omega scales, full/bulk/near views, toggles, corner-only patches,
чёрные solid walls и постоянное no-slip statement. HTML не служит источником
scientific norms. Новая metadata содержит seed ranges, S, glyph margin,
P2 wall source, оба wall maxima, profile depths/ranges/fixed scales.

## 12. What the corrected animation shows physically

**No-slip means v=0 exactly on the solid wall. It does NOT mean v=0 in a
finite-width neighborhood of the wall.** Therefore a valid viscous solution
can have v(−1,z)=0 while v(−0.99,z)≠0 and v(−0.95,z)≠0, with a large
wall-normal velocity gradient. This gradient is part of viscous
boundary-layer physics.

Практически: сначала смотрите на чёрную стенку и точку (ε=0,|v|=0), затем
на профиль внутри жидкости. В начале bulk движение у левой стороны
направлено вверх, но не означает скольжения по solid boundary. После
overshoot направление w меняется; speed не содержит знак, поэтому отдельный
diagnostic показывает **оба signed components**. Все profile/color/arrow
scales фиксированы: вязкое затухание нельзя скрыть framewise autoscaling.

Bulk-поведение не изменилось: период ≈1.61325 s; mechanical energy убывает
от 3.9843962×10⁻⁷ до 8.7189141×10⁻⁸ m⁴/s². Speed color range по-прежнему
[0,0.0011726242906210666] m/s, omega range ±0.005971933943219597 s⁻¹.
Eta quantitative panel в mm, ×1; linear tracer displacement ×181 с одинаковым
множителем x/z. Velocity решена на **fixed reference domain z≤0**; никакой
moving mesh или искусственный strip не нарисован.

## 13. What remains unresolved physically

STEP 2.1 исправляет семантику рисунка, не повышает spatial resolution и
не меняет qualification STEP 1.6. Pointwise contact slope и corner omega_max
по-прежнему не квалифицированы. Профили внутри помеченного patch показывают
дискретную модель, не точную микрофизику контакта. Не моделируются moving
contact line, wetting/dewetting, wall film, film thickness или drainage.
Нет validated nonlinear drift/mixing по траекториям линейного поля.
Fine — численный reference, не точное continuum решение.

Ни искусственная вязкость, ни slip, ни surface tension, ни smoothing/clipping,
ни подмена u/w не добавлены. Исторические STEP1_6 report/summary сохранены;
STEP2_REPORT дополнен только ссылкой на эту итерацию. Дальнейшая физика
moving-contact-line/film относится к другой модели и здесь не начиналась.

**Итог STEP 2.1: no-slip visual correction accepted.** Все шесть final MP4
прошли полное декодирование, все 130 локальных tests прошли, исходники
соответствуют green CI snapshot. Финальный provenance audit имеет
`status=passed`; оба HDF5 и исторические STEP 1.6 report/summary неизменны.
Наука ограничена прежним условно квалифицированным bulk model; контактный
угол этой визуальной коррекцией не становится разрешённой физикой.
