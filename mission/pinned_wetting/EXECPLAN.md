# EXECPLAN — pinned-wetting-PW1

## Финиш
Реальные surface data 0–5 с, необратимое связное wall wetting, final MP4,
самостоятельное уточнение и проверка результата для объявленной модели.

## Текущее состояние
**COMPLETE_FOR_DECLARED_MODEL — DECLARED_EXTENSION, 2026-09-10.**
`package --extension` session51685 завершён exit0 после112.509s. Все обязательные
gates PASS: contract, physics, wetting, numerical_refinement, rendering,
artifact_integrity и независимый review.48 итоговых artifact hashes дополнительно
проверены после упаковки;38 локальных ссылок рабочие. Активных jobs нет.

Вход: `sloshing_visualization/output/pinned_wetting/declared_extension/index.html`.
Там находятся реальные final_animation.mp4 и wall_detail.mp4 (0–5 физических секунд,
501 кадр,1920×1080,25fps), CSV histories и11 snapshots, native manifests,
FINAL_REPORT.md, verification.json и подписанный review. Исходные большие native
поля сохранены для всех пяти5s cases в results/pinned_wetting/declared_extension.

Выбран navier-space_refined-dd9d8f7c74f0:48×96,dt=.0025s,b=.50m.
P3=24.740857mm при.7875s; подъём от P1=59.661627mm.
P4=25.349900mm при2.4025s; H_L/H_R=25.349900/34.920769mm.
R_L/R_R при5s=−15.561431/+15.561431mm. Маркеры неподвижны, покрытие lower-connected
и необратимо в native/checkpoints. Fixed-b surface Linf differences:
time=7.15986e−5,space=.001390065, нормировка на начальный перепад69.841539mm.
Sensitivityb=.25/1m отдельно;60 tests PASS,0 skips. Полный независимый операторный
аудит неизменных данных выполнен реально и сохранён с source/hash; финальный
verify повторно проверил SHA и все exports/refinement/frame/full-video gates.

Это не решение исходной no-slip задачи: её solver, acceptance, reference outputs,
review и STOP сохранены. Линейная weak corner постановка не доказывает C1-гладкость
или5% физическую точность. Остаточный слой — отдельный подсеточный закон;
толщина, объём, инерция, отдельная динамика/обратное влияние не вычислялись.

Финальный ledger:4259.059063s измеренных managed jobs +1200s явно консервативного
auxiliary allowance =5459.059063s из43200s;4,907,984,082bytes новых outputs.
4GiB/job память,40GiB суммарный диск,один heavy job. Внешних расходов/install/push
не было. state.json — окончательный job-close ledger; manifest snapshot отличается
на последние0.038s закрытия job и прямо ссылается на authoritative state.

Следующее действие не требуется. Для проверки: `verify --extension` (полный
операторный проход); `resume --extension` использует завершённый пакет без PDE.
Неисполненное дальнейшее уточнение: `run --extension --mesh-level 1 --dt-scale 0.5`
создаёт отдельный output при прежних физических параметрах и ресурсных caps.

## Архив внутренних этапов до COMPLETE (не текущий статус)
Независимый physics_review завершён: PASS_DECLARED_EXTENSION,34 immutable artifact
bindings, native/audit/test SHA. Diagnostic95835 exit0,41.989261454s (в auxiliary
upper bound),8 actual MP4 frames сопоставлены с PNG,MAE<.84/255. Его jobs завершены.
Заключение extension/review.md и independent_review.json скопировано в пакет.
Следующее немедленное действие: package --extension, без PDE.

Окончательный renderer-only session90413 завершён exit0. Оба v3 MP4 полностью
декодированы,501 frame,1920×1080,25fps; кадры0/.79/2.41/5s просмотрены.
Сейчас отдельный physics_review читает окончательный пакет и готовит независимую
привязанную к SHA подпись. Его ограниченный read-only diagnostic session95835:
пять native SHA, eta/H/CSV/refinement и8 декодированных MP4 кадров vs PNG;
без нового PDE/FEM/operator solve, учёт в auxiliary allowance. Активных heavy jobs нет.
После завершения этого diagnostic и подписи автоматически
`package --extension` проверит native hashes и все export/refinement/render gates.
Перед этим budget:4146.550s измеренных jobs +1200s auxiliary upper bound =5346.550s
из43200s; новых файлов4,907,955,319bytes. Исходные no-slip milestones восстановлены
в поле reference_milestones из неизменённого исторического STOP.

Все пять DECLARED_EXTENSION расчётов 0–5 s завершены. Полный независимый аудит
каждого native velocity/pressure/RK stage завершён: scientific_checks_passed=true,
пять native_reports PASS. Session60525 завершён exit2 только из-за ожидаемой
отдельной review-подписи; kernel job lock свободен, active_job=null.
Выбранный run: navier-space_refined-dd9d8f7c74f0,48×96,dt=.0025s,b=.5m.
P3=24.740857mm при.7875s; P4=25.349900mm при2.4025s; H_R=34.920769mm.
Следом: общий60-test suite; renderer-only повтор для читаемых подписей близких
P3/P4 и фиксированного увеличенного wall panel; независимый visual/native review;
package --extension. Уже выполненный полный native audit фиксируется с исходником
и SHA всех данных; повторное использование требует отдельной reviewer-подписи.
Все export/refinement/film/frame/full-decode проверки при упаковке идут заново.
Обычный CLI verify --extension по-прежнему пересчитывает полный native audit.
Общий suite завершён:60 PASS,0 skipped,27.71s (extension/tests.xml).
Renderer-only session68856 завершён exit0; v2 кадры просмотрены мной и reviewer.
Перед последним рендером добавлены требуемые RENDERING.md горизонтальные пунктирные
ориентиры исходных P1/P2 и явная подпись разных масштабов x/z. Координаты неизменны;
v2 автоматически архивируется. Физические source/native файлы не менялись.
Учтён дополнительный600s auxiliary
allowance, суммарно1200s отдельно от измеренных heavy job интервалов.
Ограничение corner: движущийся endpoint не совместим с C1-напряжениями ровно
в пересечении free-stress/Robin граней; weak BC на открытых гранях допускают
угловую нерегулярность. Uniform C1 convergence и 5% физическая точность не заявлены.

## Журнал предыдущих внутренних этапов
После исходного STOP пользователь явно разрешил отдельный DECLARED_EXTENSION:
Navier slip боков, непроницаемость, no-slip дна, sigma=0 и независимый retained layer.
Прежний STOP архивирован в history/no_slip_STOP_state.json; его данные и acceptance
не меняются. Новый model_decision/CONFIG/ACCEPTANCE_ADDENDUM — в extension/.
b=0.50 м выбрана до pilot, чувствительность 0.25/1.00 м. Изготавливается целевой
путь с реальными свободными вертикальными endpoint DOFs и положительным трением.

Pilot `navier-pilot-3ef311facf94` принят: 201 state, 0–1 s, max eta/a=0.034921,
max exact wall-inclusive slope=0.082026. Channel benchmark прошёл. Preview MP4
101 frame, 1920×1080,25fps полностью декодирован; видимая полоса согласована с H.
Текущий следующий job: baseline24×48 dt=.005 5s, time24×48 dt=.0025,
space48×96 dt=.0025, затем b=.25/1 на fine mesh. Один kernel-locked job.
Для accepted RK evidence лимит одного файла увеличен с 1 до 8 GiB;
общий неизменный лимит 40 GiB, RAM 4 GiB/job и 12h остаются обязательными.
Независимый verifier реализуется отдельным reviewer context без импорта Navier solver.

2026-09-10: прочитаны все нормативные документы. Ветка уже
`feature/pinned-wetting-history`, HEAD совпадает с observed master `7f0c45f`.
Пользовательские untracked AGENTS.md и mission/ сохраняются. Реальный reference
`reference-20260909T222724-02a9a3` завершён до 1 с: 401 accepted state, preview MP4
101 кадр, 1920×1080, 25 fps, полное декодирование прошло. Целевой расчёт NOT RUN.
Независимый review выявил недостающее физическое замыкание и помог усилить verifier.

## Milestones
- [x] M0: audit + resolved_case/model_decision + сквозной reference preview; target model blocked.
- [x] M1 расширения: свободный wall endpoint с конечным Navier slip и отдельный retained state; pilot 1s и preview готовы.
- [x] M2: target 5 s + предварительное видео + отдельное уточнение dt и mesh + b sensitivity.
- [x] M3: окончательные MP4/data/report + independent verify; package exit0, COMPLETE_FOR_DECLARED_MODEL.

Baseline расширения `navier-baseline-6b042d9d8353`: 1001 states до5s; локальные
рекорды P3=24.746841mm при.785s, P4=25.370608mm при2.400s. Это ещё предварительные
значения до независимого mesh/dt сравнения, не storyboard. 6 ограниченных Navier
тестов PASS: нулевое равновесие, три b channel benchmarks, реальные свободные endpoints,
bitwise checkpoint continuation после рождения P3 (включая все RK поля и coating).

Time refinement `navier-time_refined-03e9ce1500fc` завершён 0–5s/2001states.
При фиксированных b=.5 и24×48: normalized surface Linf=7.15986e-5,
R Linf=7.06022e-5,H Linf=7.82852e-6; first-peak time difference fraction=.00317460.
Первый peak=.7875s,24.747252mm; P4=2.4025s,25.371155mm. Все time-gates пройдены.
Текущий heavy job60525/PID2 в своём namespace выполняет space48×96; state.json
хранит актуальное физическое время и kernel-lock identity. Нельзя запускать дубликат.
Общий ограниченный pytest:46 PASS,0 skipped; tests.xml в extension/. Два предупреждения
только из грубого historical reference continuation test. Новые Navier проверки PASS.
Визуально обнаружено наложение подписей близких P3/P4 в preliminary; запланирован
один адресный renderer-only повтор с text offsets/выносками. World coordinates не меняются.

Space `navier-space_refined-dd9d8f7c74f0` завершён:2001states0–5s.
P3=24.740857mm,.7875s;P4=25.349900mm,2.4025s;оба уточнения прошли gates.
Max slope вырос.082026→.103589, convective indicator.110717→.125357;
не заявлять сходимость угловых производных или5% физическую точность.
Wall traction L2 residual уменьшился1.30724e-4→5.58140e-5. Линейные screens пройдены.
Следующие автоматические cases — b=.25/1m на48×96 dt=.0025, без изменения выбранногоb=.5.

55 тестов PASS,0 skips (extension/tests.xml). Дополнительный allowance600s для
коротких extension tests/diagnostics будет списан при следующем CLI job, сверх
сохранённых600s reference allowance; в extension.py есть одноразовый флаг.
Это консервативный upper bound, не измерение. Все heavy solver/render/verify jobs
учитываются измеренно. Ни один фактический cap не близок к исчерпанию.

b=.25 run `navier-slip_low-62f92521a3d1` завершён2001states0–5s;
P3=24.066818mm,.795s;P4=24.387649mm,2.410s.
Сейчас тот же job60525 выполняет b=1 `navier-slip_high-6eaa7fb099eb`.
После него: final render -> independent verifier в текущем job. Detector sensitivity
и record_stability уже созданы из неизменённых трёх fixed-b datasets; их повторный
postprocess не меняет PDE или H. Исходник renderer до адресного исправления
сохранён в extension/history/renderer_before_label_fix.py (SHA совпадает).

Все пять5s cases завершены. b=1 run `navier-slip_high-6eaa7fb099eb`:
P3=25.114557mm,.7825s;P4=25.888842mm,2.3975s. Max surface differences vsb=.5:
b=.25→.0191855, b=1→.0107486 (нормировка на начальный перепад, не физическая ошибка).
Все amplitude/slope screens пройдены. b=.25: max slope=.179689, convection/A0=.184360,
kinematic/U0=.118671 — ограничения заметнее, чем у выбранногоb=.5.
Job60525 сейчас RENDERING, затем проверит native/exports независимым verifier.
Не менять extension_output.py до завершения текущего рендера: процесс держит старый
модуль в памяти; signature должна соответствовать фактически использованному коду.
После статуса VERIFYING можно править только renderer; после освобождения job.lock
выполнить `render --extension` (без PDE), затем вызвать physics_review для итоговой
raw/visual подписи и `package --extension`. На момент этой записи renderer-only
исправление ещё НЕ выполнено. Все numerical source files остаются неизменными.

## Историческое следующее действие (до выпуска)
Текущее действие расширения: `python sloshing_visualization/scripts/pinned_wetting_mission.py run --extension`.
Внутри: pilot -> target5s -> отдельные dt/mesh -> slip sensitivity -> final render -> независимый verify.
Историческое описание предыдущего STOP ниже сохраняется как история, не текущий blocker.
Дальнейший target требует одного физического решения из model_decision.md.
Все доступные независимые части завершены, оба reference run сохранены.
Последний `reference-20260909T223914-831392` прошёл stop-after40/resume до400,
битwise совпал по физическим полям с первым; source snapshots раздельные.
Итоговый вход: `sloshing_visualization/output/pinned_wetting/index.html`.
Независимое заключение: review.md / independent_review.json,
`PASS_REFERENCE_REJECT_TARGET`. Target 5 s/refinement/final MP4 остаются NOT RUN.

## Реестр решений
2026-09-10: выбраны лимиты CONFIG: 43200 s совокупного локального job wall time,
40 GiB новых файлов, не более 65% host RAM, один тяжёлый процесс.
Host: 8 logical CPUs, 33.25 GB RAM, около 44.57 GB свободного диска.
Runtime переиспользован read-only из
`/home/victor/Projects/waves-analyzer/sloshing_visualization/.venv/bin/python`.
Новая установка не выполнена: sandbox network failure, попытка approval прервана,
после чего найдена уже установленная dependency в другом worktree.
Бюджет jobs ведётся фактическими интервалами в state.json, lock ядра исключает
дубли даже между разными PID namespaces. RAM дополнительно ограничена 4 GiB/job.
До итоговой упаковки/проверок управляемые jobs заняли 80.106 с. Для коротких
диагностических и тестовых команд отдельно консервативно списано 600 с;
это верхняя бюджетная оценка, не притворное измерение. Итого списано 680.106/43200 с.
Новые файлы около 67 MiB из 40 GiB; без внешних расходов и push.
Reference использует неизменённый linear P2/P1 FEM; его fixed endpoints не выдаются
за moving contact. Данные и видео получат REFERENCE / NOT TARGET.
Остаточный слой: аналитическое связное множество [bottom,H] с нулевым объёмом
в ведущем порядке; толщина/инерция/обратное влияние/adsorption energy не вычисляются.

## Учёт jobs
Авторитетные PID/start ticks/namespace, kernel lock и heartbeat находятся в state.json.
Не доверять устаревшей записи плана без проверки kernel lock.

## Выполненные команды / evidence
Прочитаны FEM trace R, integrators и Diagnostics; endpoint rows R тождественно нулевые.
Прочитан require_scope full_phase_rate: matched density, zero force; guard сохраняется.
В текущем worktree исторических HDF5 нет; найдены в соседнем waves-analyzer,
но не импортированы как новая физическая траектория. Исторические данные не изменены.
Независимый reviewer: /root/physics_review (read-only физический аудит, затем data review).
46 тестов PASS, включая 12 обязательных негативных категорий и 10 мутаций реальных
экспортированных файлов с обновлёнными checksums. Evidence: test_results.xml.
Точное bitwise продолжение checkpoint, нулевое hydrostatic disturbance, lock,
coating connectivity и неподвижные marker IDs проверены.
Reference-аудит: weak divergence 1.87518e-14, относительная масса 1.27719e-17,
wall velocity 0, энергетическая невязка 2.87607e-14. Это нормы reference, не target PASS.
Максимальный локальный наклон reference 9.33957: явно непригоден для квалификации
пристеночного контакта линейной геометрией.
`package` выпустил верхнеуровневые manifest/verification/FINAL_REPORT и независимый
review без изменения закрытых native/run bundles; return code 2 означает model blocker.
`run` с прежними контролами и `resume` повторно используют готовый reference
после read-only проверки; лишний PDE не запускается.
Последний контроль `run` вернул `PDE_restarted=false`; список run IDs остался из двух.
20 локальных HTML/Markdown ссылок проверены, битых нет; handoff hashes PASS.
`verify --reference-only` exit 0; обычный `verify` exit 2, target acceptance не пройден.
Финально учтено 82.942 с управляемых jobs, с отдельным консервативным auxiliary charge
682.942 с из 43200 с. Новые данные 70,246,082 bytes. Активных jobs нет.

## Блокеры
Следующий абзац — сохранённый исторический blocker исходной no-slip задачи,
не причина остановки разрешённого пользователем DECLARED_EXTENSION.
Один blocker PW1_CONTACT_CLOSURE: нужен закон нового контакта и его масштаб,
согласованные с no-slip, sigma=0, массой и энергией. Не нужна толщина слоя.
Резкая материальная поверхность даёт нулевые endpoint rows и нулевое движение
концов; sigma=0 обнуляет имеющиеся CHNS interfacial/wall energy terms.
Диффузный no-slip контакт в принципе возможен при заданном замыкании.
Команда: `python sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check`.
Подробности: model_decision.md; два пути проверены без изменения historical guards.

## Перед сокращением контекста/завершением сессии
Обновить state.json, реальные job IDs, текущий checkpoint, бюджет и exact next command.
