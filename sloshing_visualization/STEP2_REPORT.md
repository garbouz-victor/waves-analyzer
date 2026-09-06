# STEP 2 — scientific bulk-flow visualization

## 1. Goal

**CONDITIONALLY QUALIFIED FOR BULK MODEL VISUALIZATION.**
Показать рассчитанную динамику bulk жидкости, обмен энергии, смену направления
скорости и интегрально согласованные вихревые структуры. Это новый слой
интерпретации, не изменение глобального verdict STEP 1.6. Никакая новая
production PDE simulation в STEP 2 не запускалась. Solver, FEM, физические
граничные условия и исторические STEP 1.6 report/summary не изменены.

## 2. Source dataset

Источник численного решения — commit
`dc51e33fc18f598e0db311d68d6cdcce3d9b68ac`, «step 1.6 qualify physical animation dataset».
Основной HDF5: `validation_results/animation_qualification/runs/fine-full.h5`.
Medium из того же каталога — только comparison/provenance, не rendering source.
Оба прошли существующий `validate_cache()` и переиспользованы без изменения.

| Параметр | Значение |
|---|---:|
| a, d, g | 1 m, 10 m, 9.81 m/s² |
| alpha, nu | 0.02°, 0.01 m²/s |
| fine intervals / triangles | 120×210 / 50400 |
| medium intervals / triangles | 80×140 / 22400 |
| x/z grading | 1.8 / 4.0 |
| velocity / pressure / eta / omega | P2 / P1 / P2 trace / DG1 |
| integrator / dt / snapshot_dt | SDIRK2 / 0.00125 s / 0.005 s |
| physical time / snapshots | 0…5 s / 1001 |

`nx=null,nz=null` в исходном config означает выбор preset, а не неизвестную
сетку: фактические интервалы независимо подсчитаны из mesh coordinates.

| HDF5 | Bytes | mtime_ns |
|---|---:|---:|
| fine | 2364886404 | 1788665248766948120 |
| medium | 1072567328 | 1788657587171233015 |

Полные SHA-256 содержимого (не путать с hash только config):

```text
fine   fa299d74c9487a40a5fac55966c1288c950bc0a46116f80e7da0063f794cf09a
medium 1860e6e07236e2bab9511fd0aea876ca7e57ddf5cfee019522499fc1a08bc94a
```

Revision fingerprint для кэширования включает path, bytes, mtime_ns и
config SHA-256. Все эти поля, полный content hash, resolved mesh и хэши
кода renderer сохранены в [artifact_audit.json](output/step2/artifact_audit.json).
Факторы, сетки display, версии библиотек и времена событий — в
[animation_metadata.json](output/step2/animation_metadata.json).

## 3. What STEP 1.6 established

Исторический результат **NOT QUALIFIED FOR PHYSICAL ANIMATION** сохраняется:
глобальный slope fine=0.373408, первое сохранённое превышение 0.3 — t=2.365 s.
Это локальная contact structure, не необъяснённая неустойчивость solver.

| Измерение STEP 1.6, 0…5 s | Результат |
|---|---:|
| medium/fine bulk eta difference | 0.075961% |
| medium/fine velocity difference | 0.079293% |
| wall omega field difference | 5.73043% |
| bulk path mesh separation | 0.402706 µm |
| no-slip / contact residual | 0 / 0 |
| fine max volume error | 1.103×10⁻¹⁹ m² |
| fine relative energy budget residual | 4.520×10⁻¹³ |
| max relative strong divergence, medium / fine | 8.428% / 6.443% |
| fine bulk max slope, abs(x)≤0.90 | 4.53719×10⁻⁴ |
| fine intermediate max slope, abs(x)≤0.98 | 4.08995×10⁻³ |

Field differences здесь нормированы на peak fine norm за весь интервал,
не на почти нулевое поле отдельного кадра. Это не новые нормы display grid.
Сильная дивергенция не исчезает от слабого Bv≈0: Taylor–Hood не является
pointwise divergence-free. Fine уменьшает её, но остаётся дискретным приближением.

SHA-256 исторических файлов до/после STEP 2 совпадает:

```text
STEP1_6_REPORT.md
85da88bb6726ba486c2048d9a9048d6bfc88e65b0a9ee672f20ff3c0fcd456ef
validation_results/animation_qualification/summary.json
997266bcecc6a2a21947119dba5c9d132571b423d27e832547365b91e4157858
```

## 4. Why no new smaller-alpha run is needed for bulk visualization

Production PDE линейна; STEP 1.6 независимо проверил масштабирование по
tan(alpha). Уменьшение alpha уменьшает абсолютное смещение, но не relative
medium/fine errors, relative divergence или локализацию pinned-contact
компонента при сгущении. Новый alpha=0.0025° run не сделал бы continuum
corner разрешённым. Здесь сохранён реально рассчитанный alpha=0.02° dataset;
не применяется даже amplitude rescaling численных полей.

## 5. Fixed-domain interpretation

Velocity is solved on the fixed reference domain z≤0. The plotted free
surface is the first-order interface displacement. Mesh не движется.
Поле, стрелки и их наконечники остаются в reference domain; поверхность
рисуется отдельной линией. Скорость не продолжена в z>0 и не деформирована
вместе с eta. Нижняя граница near-surface окна −1.5 m — граница zoom,
а не дно сосуда (дно находится на −10 m).

## 6. Qualified bulk region

Machine-readable convention: [animation_scope.json](configs/animation_scope.json).
Bulk |x|≤0.90; transition/caution 0.90<|x|≤0.98. Это визуальные зоны,
основанные на STEP 1.6, не математическая граница применимости в x=0.98.
Основная интерпретация: bulk surface/velocity, reversal, энергия, затухание,
penetration и tracers вдали от углов. Transition слегка оттенён.

## 7. Unqualified contact region

0.98<|x|≤1 показана серой штриховкой в explained/clean/vorticity/surface/
tracer videos. Контакты отмечены неподвижными маркерами. Постоянная подпись:
«Pinned-contact region: not physically resolved by this model».
Eta и omega здесь не сглажены, не обрезаны и не заменены плёнкой. В HTML
можно выключить заливку, но dashed границы и предупреждение остаются.
Рост corner slope/omega_max при сгущении не выдаётся за точную контактную физику.

## 8. Display magnification policy

| Величина | Политика, неизменная на всех кадрах |
|---|---|
| eta в flow panel | ×1, metres |
| eta в количественной панели | mm, только перевод единиц; пределы ±0.390954 mm |
| tracer glyph | Xdisplay=X0+181·(Xlinear−X0), одинаково x/z |
| speed background | physical m/s, без amplification |
| arrows | fixed factor 102.334568 s; key 0.5 mm/s → 0.0511673 m рисунка |
| geometry flow/tracer panels | equal aspect; никаких разных x/z magnifications |

True eta extrema за dataset: ±0.349065865 mm. Они не затухают как global
extrema из-за pinned endpoints; затухание видно в bulk профиле и modal signal.
Максимум физического linear tracer displacement=0.527693 mm, p95=0.380980 mm.
Первоначальный target p95 display=0.075 m ограничен одним глобальным фактором
по запасу до границ для **обеих** tracer models; после округления M=181.
Итоговый p95 glyph displacement≈0.06896 m. Никакого покадрового изменения M,
анизотропии, отражения или clipping траекторий нет. Pure transforms не меняют raw arrays.

## 9. Velocity rendering

P2 FEM maps построены один раз для фиксированной near grid 161×121 и full
grid 81×161. Каждый snapshot читается отдельно; нет загрузки 2.36 GB в RAM.
Display cache около 295 MiB, velocity/omega float32; eta/tracers/histories
float64. Сохранённые coefficients не изменены. Повторная FEM оценка семи
key/end snapshots дала максимальные разности округления u=2.91×10⁻¹¹ m/s,
w=5.76×10⁻¹¹ m/s, omega=9.91×10⁻⁹ s⁻¹; eta совпала точно.

Speed range **[0, 0.0011726242906210666] m/s** — максимум fixed near grid
в |x|≤0.98 по **всем** 1001 snapshots, а не максимум текущего кадра.
Это display statistic, не новый FEM error benchmark. Quiver 25×18,
линейная постоянная шкала; нулевая скорость не превращается в minimum-length dots.
Exact-wall sampling на каждом snapshot дал **0 m/s**. Все arrow tips:
−1≤x≤1, max z=−0.0544366 m. На heatmap заданы явные half-cell boundaries:
верхний pixel не выступает выше z=0.

## 10. Surface rendering

Профиль вычисляется из P2 trace на каждом surface edge, включая graded
contact edges: 961 drawing points, не грубый uniform spline. Глобальный
диапазон дополнительно проверяется по точным стационарным точкам P2.
Показаны eta=0, initial dashed line, текущий профиль и pinned endpoints.
Surface panel всегда количественная, в mm. Её разный масштаб осей x и eta
подписан единицами; это график, не геометрия увеличенного сосуда.

## 11. Energy interpretation

K/P/E берутся непосредственно из HDF5 diagnostics и синхронизируются по
snapshot time. Единицы: m⁴/s², на единицу плотности и поперечного размаха.
E(0)=3.984396174×10⁻⁷, E(5)=8.718914113×10⁻⁸ m⁴/s².
Необъяснённого роста нет. Вязкая диссипация рассчитана PDE; RK correction
учтена в исходном энергетическом бюджете, **не** изображается физическим теплом.
Первый/global K peak t=0.395 s, K=3.474744489×10⁻⁷; это близко, но не равно
first modal crossing 0.416098 s. Ни моменты, ни равенство фаз не подгонялись.

## 12. Vorticity rendering and its limits

Omega=dw/dx−du/dz интерполируется из сохранённого **DG1 FEM field**, не из
finite differences display velocity. Цветовая шкала симметрична и постоянна:
**[−0.005971933943219597, +0.005971933943219597] s⁻¹**.
Это 99.5 percentile |omega| в |x|≤0.98, every second near-grid row/eligible
column на всех snapshots. Пределы вычислены однократно. Вне диапазона
значения **сохраняются**, цвет насыщается с extend indicators на colorbar.
Исключение corner extrema из color-scale statistics не удаляет их из поля.

Integral wall-vorticity structure numerically coherent на измеренном
~5.7% medium/fine field-difference уровне. Видны слои у no-slip стенок,
смена знака и перестройка при reversal. Pointwise corner extrema не разрешены;
пиксели и display percentile не измеряют физический размер контактного вихря.

## 13. Linear Lagrangian tracers

Основной режим: ξ(X0,t)=∫v(X0,τ)dτ, Xlinear=X0+ξ, где **X0 фиксирована**.
Скорость берётся непосредственно P2 FEM evaluation в seed points.
Трапеции дают точный интеграл выбранного piecewise-linear temporal interpolant
и второй порядок для smooth synthetic signal, проверенный тестом.
Никакого Euler stepping нет. Main показывает 21 tracer на z0=−0.1,−0.3,−1;
comparison показывает все 28, включая z0=−3. Цвет задаёт исходную глубину;
крест — X0, tail — последние 0.5 s.
Проверка reference domain и reconstructed true surface не обнаружила invalid
trajectories у выбранных seeds. Invalid не отражаются и не clip'аются.

## 14. Difference between linear tracers and advected pathlines

В comparison справа решается dX/dt=vlinear(X,t) через прежний RK4 + P2 FEM
и интерполяцию между snapshots. Слева — first-order displacement.
Обе панели имеют одинаковый M=181. Максимальная истинная разность
**0.313114 µm** за 0…5 s, на t=5 — около 0.2975 µm.
Spatially uniform field даёт совпадение (тест включает настоящий RK4 path);
при varying field возникает higher-order разность. Оценка скорости уже в
смещённой X вводит O(alpha²) компоненты без отсутствующей Eulerian v2.
Ни net drift, ни mixing из этих paths не считаются физически проверенными.

## 15. Depth penetration

Синхронный full-depth профиль Q(z)=sqrt(∫|v(x,z)|²dx), −10≤z≤0, получен
из горизонтальной FEM quadrature. Q имеет единицы m^(3/2)/s.
В panel depth — горизонтальная ось, Q — подписанная log-ось с фиксированными
пределами. Нулевые значения не подменяются положительным floor (не видны на log).
Это профиль, а не сжатое изображение сосуда.

В t=0.395 s отношения Q(z)/Q(0):

| фактический z, m | отношение |
|---:|---:|
| −1.0125 | 0.203521 |
| −2.916667 | 0.00855518 |
| −5.041667 | 0.000247052 |
| −10 | 0 |

Так глубина 10 m остаётся видна, хотя главное окно фокусируется на верхних 1.5 m.
Explorer дополнительно предоставляет full-depth heatmap с equal aspect.

## 16. Modal period and phase markers

Modal projection sin(πx/(2a)) — диагностический функционал рассчитанной eta,
не потенциальная модель вместо PDE. События заново получены прежней проверенной
функцией `modal_events` из fine diagnostics; renderer не содержит заданных времён.

| Событие | Измерено, s | Ближайший snapshot, s |
|---|---:|---:|
| release | 0 | 0 |
| first zero | 0.416098 | 0.415 |
| opposite extremum | 0.809794 | 0.810 |
| next zero | 1.222442 | 1.220 |
| next positive extremum | 1.615950 | 1.615 |
| late opposite extremum | 4.037690 | 4.040 |

Измеренный characteristic period **1.61324728 s**, same-sign amplitude ratios
0.77756…0.77992 за цикл. Inviscid reference 1.60060963 s служит только сравнением,
не генератором кадров. Phase labels кратко описывают измеренное событие/интервал.

Знаки подтверждаются непосредственно calculated samples. В t=0.1:
u(0,−0.1)=−0.301764 mm/s, w(−0.75,−0.1)=+0.343268 mm/s,
w(+0.75,−0.1)=−0.343268 mm/s. В t=1.0 знаки меняются:
u(0,−0.1)=+0.646278 mm/s, w слева/справа=−/+0.474162 mm/s.
Последовательность release → переток → crossing → overshoot → reversal →
затухающие циклы происходит из полей; симметрия и энергия не дорисованы.

## 17. Animation outputs

Все пять MP4: **1920×1080, 50 fps, 1001 кадр, 20.02 s**. Использованы все
сохранённые состояния без промежуточных synthetic PDE frames. Замедление ×4;
последний snapshot удерживается один frame interval 0.02 s.

| Артефакт | Размер MB (decimal) | Encoding runtime, s |
|---|---:|---:|
| [bulk_flow_explained.mp4](output/step2/bulk_flow_explained.mp4) | 5.109 | 445.7 |
| [bulk_flow_clean.mp4](output/step2/bulk_flow_clean.mp4) | 4.941 | 459.7 |
| [bulk_vorticity.mp4](output/step2/bulk_vorticity.mp4) | 3.699 | 448.2 |
| [free_surface_true_scale.mp4](output/step2/free_surface_true_scale.mp4) | 1.205 | 165.9 |
| [tracer_model_comparison.mp4](output/step2/tracer_model_comparison.mp4) | 0.934 | 255.1 |

Дополнения:

- [key_phases.pdf](output/step2/key_phases.pdf): 9 страниц, события и поздние
  extrema, actual time и численные K/P/E; [PNG contact sheet](output/step2/key_phases.png).
- [explorer.html](output/step2/explorer.html): автономный Plotly, 101 state,
  near 81×61/full 41×81, 36.6 MB, генерация 11.6 s. Slider/play/pause, поля,
  стрелки, tracers, mask, виды; шесть significant digits для display fields.
  Нет CDN; browser timing приблизительное, не заменяет cadence MP4.
- [t0](output/step2/frames_preview/t0.png),
  [first zero](output/step2/frames_preview/first_zero.png),
  [opposite extremum](output/step2/frames_preview/opposite_extremum.png),
  [next zero](output/step2/frames_preview/next_zero.png),
  [next positive](output/step2/frames_preview/next_positive.png),
  [late cycle](output/step2/frames_preview/late_cycle.png).
- [preview_metrics.json](output/step2/preview_metrics.json),
  [render_summary.json](output/step2/render_summary.json),
  [inspection_summary.json](output/step2/inspection_summary.json),
  [explorer_inspection.json](output/step2/explorer_inspection.json),
  [test_summary.json](output/step2/test_summary.json).

Static preview gate: шесть обязательных кадров + четыре варианта просмотрены
до полного rendering. В таблице скорости — near display grid, не новый FEM norm;
энергии — реальные HDF5 diagnostics в единицах 10⁻⁷ m⁴/s².

| t, s | speed max, mm/s | bulk speed max, mm/s | K | P | E |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 3.984396 | 3.984396 |
| 0.415 | 1.124788 | 1.117752 | 3.457226 | 0.058740 | 3.515966 |
| 0.810 | 0.371290 | 0.239311 | 0.037583 | 3.094878 | 3.132461 |
| 1.220 | 0.948269 | 0.948269 | 2.720485 | 0.039535 | 2.760020 |
| 1.615 | 0.435301 | 0.276208 | 0.037863 | 2.406198 | 2.444061 |
| 4.040 | 0.277944 | 0.168252 | 0.012972 | 1.138101 | 1.151074 |

Во всех этих кадрах exact eta range ±0.349066 mm (pinned endpoints), speed
minimum=0, M_eta=1 и M_tracer=181. Дополнительный short preview: 0…1 s,
201 frame, 4.02 s, 1.171 MB, encoding 78.8 s; он просмотрен перед full outputs.

### Validation of the visualization itself

**117 ordinary tests passed (16.65 s), 2 validation tests passed (26.94 s).**
Итого 119; все прежние STEP 1/1.5/1.6 тесты сохранены. 23 новых STEP 2 tests
проверяют scope/events, source maps, exact walls, energy field alias,
fixed color limits на разных амплитудах, drawing boundaries, contact mask,
pure magnification, linear tracer order/zero/uniform/higher-order cases,
fingerprint invalidation, preview gates и explorer contract.
Полный перечень — в test_summary.json. Обычный suite вывел 15 warnings:
14 ожидаемых diagnostics на coarse validation cases и один fontTools deprecation.
Validation suite: один fontTools deprecation. Warnings не подавлялись.

Все пять полных видео и short preview полностью декодированы ffmpeg;
ffprobe подтвердил кадры, размер и каждый PTS interval=0.020 s. Просмотрены
извлечённые key/end кадры, storyboard и реальный browser screenshot.
Browser smoke test проверяет slider, play/pause, fixed speed/omega limits,
три view, toggles и постоянное contact warning; runtime exceptions отсутствуют.
Это проверка корректности артефактов, не новая physics qualification corner.

Исправленные проблемы разработки не скрыты:

1. Ошибочное имя potential-energy reader заменено на исходное
   `gravitational_potential_energy`; добавлен regression test.
2. Depth inset закрывал часть стенки: перенесён в отдельную правую панель.
3. Matplotlib nearest shading продолжал верхнюю половину pixel в z>0.
   Первый full render остановлен; введены explicit cell edges, тест;
   static/short gates повторены, затем все full outputs пересозданы.
4. Plotly сохранял aspect constraint старого view: view switch теперь явно
   пересоздаёт layout; browser regression подтвердил правильные ranges.

Неудачный initial display cache и superseded preview/partial video сохранены
локально в ignored quarantine paths, не переиспользованы. Источники не тронуты.

### Reproduction and Git hygiene

Точные команды и preview gates приведены в [README](README.md#step-2--scientific-bulk-visualization).
Renderer не запускает solver и отвергает missing/unfinished/incompatible source
или cache. При изменении source revision нужно явно создать новый display output.
FFmpeg через [Matplotlib FFMpegWriter](https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FFMpegWriter.html);
HTML использует [Plotly heatmap с явными edges](https://plotly.com/javascript/reference/heatmap/).
Precompute занял 122.3 s; это FEM postprocessing, не PDE run.

В Git: source, config, tests, metadata/audit, report, 10 preview PNG и key PDF/PNG.
Не в Git: HDF5, display cache, MP4, HTML, decoded audit frames, raw pytest XML.
Полные MP4 занимают суммарно 15.89 MB; несмотря на умеренный размер, они
regenerated artifacts и оставлены ignored для чистой истории. После clone
локальные video/HTML ссылки требуют воспроизведения. Битовая идентичность MP4
на иной версии ffmpeg не обещается; physics provenance и display policy явны.

## 18. How to watch the animation

1. Сначала смотрите на TRUE-SCALE eta panel в mm: bulk наклон, crossing,
   overshoot; endpoints остаются неподвижными.
2. Затем одновременно смотрите на velocity arrows: после release вверху
   течение справа налево, слева подъём, справа опускание.
3. Сравните с kinetic-energy curve: K достигает первого максимума около
   прохождения modal disturbance через equilibrium, не ровно в тот же момент.
4. Следите за tracer layers ×181 и Q(z): движение резко уменьшается с глубиной.
5. После overshoot смотрите на reversal arrows, затем затухающие циклы.
6. В vorticity video следите за рождением и сменой знака wall vorticity.
7. Серые полосы у x=±1 **не** интерпретируйте как физически разрешённую
   contact-line dynamics. Даже большая насыщенная corner pixel не является доказательством.

## 19. What can be inferred physically

В рамках текущей линейной pinned-contact модели можно интерпретировать
bulk free surface, направления скорости/циркуляции, обмен K/P, вязкое
затухание, проникновение с глубиной, first-order displacement вдали от угла
и интегральную/качественную no-slip vorticity structure. Начальное поле
скорости нулевое; движение создаёт gravity/stress imbalance, а не сценарий renderer.
Модель и измеренная сеточная согласованность остаются условиями этих выводов.

## 20. What cannot be inferred physically

Нельзя считать converged corner slope, физически точным corner omega_max,
реальной толщиной плёнки локальную контактную структуру, или fine solution
точным continuum limit. Реальное закрепление контакта в эксперименте не
следует из математического no-slip выбора. Не смоделированы moving contact
line, wetting/dewetting, film deposition, drainage или нелинейное перемешивание.
Sigma=0 в PDE не доказывает малость surface tension в реальной microphysics.
Accumulated O(alpha²) pathline drift не прошёл physical validation.

Для будущего **STEP 3**, если цель — уход жидкости с высокой стенки, плёнка
и её стекание, нужна отдельная физическая ветка: moving-contact-line
regularization (обоснованный slip, precursor film либо diffuse interface),
wetting/contact-angle parameters, surface tension с проверкой масштаба,
интерфейс, допускающий не-single-valued форму/вертикальную плёнку, и nonlinear
free-boundary treatment. Потребуются новые benchmarks и независимая проверка
пространственного/временного разрешения плёнки. Ничего из этого в STEP 2 не добавлено.
