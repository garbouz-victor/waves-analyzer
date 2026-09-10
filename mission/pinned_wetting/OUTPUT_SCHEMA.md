# Конечные данные и независимая проверка

Выходная папка: `sloshing_visualization/output/pinned_wetting/`.
Native field data могут лежать в `results/pinned_wetting/<run_id>/` и быть указаны
в manifest. Все пути в итоговом пакете должны реально разрешаться; не ссылаться
на sandbox или файлы, существующие только в другом рабочем окружении.

## Минимальная структура

```
index.html
final_animation.mp4
wall_detail.mp4
resolved_case.json
model_decision.md
surface_history.csv
wetting_history.csv
marker_catalog.csv
marker_tracks.csv
runup_summary.json
energy_history.csv
refinement.json
key_frames.png
FINAL_REPORT.md
manifest.json
verification.json
```

`index.html` — рабочий вход в пакет. Можно встроить видео и ссылки на CSV, PNG,
отчёт; сложный explorer опционален. Не нужен сервер для базового просмотра.

## surface_history.csv

```
time_s,point_index,x_m,z_m,component_id,state_id
```

Кривая может быть параметрической — не навязывать z=eta(x), если поверхность
становится вертикальной. component_id различает отдельные контуры при их наличии.
Для t=0 точки принадлежат заданной прямой. Отдельно фиксируется связь области,
контуров и стенок. Редкие frame snapshots не используются для интегральной памяти.

## wetting_history.csv

```
time_s,R_L_m,R_R_m,reach_L_m,reach_R_m,H_L_m,H_R_m,state_id
```

R — мгновенная макроскопическая высота на конце принятого интервала.
reach — максимум достигнутого на этом интервале, включая проверенные внутренние
стадии/события. При отсутствии внутреннего peak finder reach=R.
Поэтому reach>=R, а H_i[n]=max(H_i[n-1], reach_i[n]).
В начальной строке t=0: R=reach=H=z1/z2.

Таблица строго возрастающая по времени, все поля конечны. state_id указывает
на принятый run/state; для reach из внутренних точек нужна отдельная source-map
в manifest. Не переносить rejected stage peak в H до принятия шага.

W_i=[z_bottom,H_i] — подсеточная связная покрытая область. Проверка таблицы доказывает
только логику памяти, не достоверность R/reach; их нужно независимо пересчитать
из native data/контактного алгоритма. Bundled checker прямо возвращает
`scientific_status: NOT_VALIDATED_BY_THIS_CHECK`.

## marker_catalog.csv

```
marker_id,side,height_m,created_at_s,state_id
```

P1: side=L, height=z1, created=0.
P2: side=R, height=z2, created=0.
Новые marker IDs создаются только для разрешённых новых рекордных локальных пиков.
Они не обязаны существовать. Нельзя переиспользовать ID для новой высоты.

## marker_tracks.csv

```
time_s,marker_id,side,height_m
```

Таблица соответствует отрисованным frames. На каждом frame присутствуют все уже
созданные маркеры, даже если они ниже текущего H. Координата/сторона ID неизменны.
Здесь нужны data coordinates, а не лишь pixels. Дополнительно renderer может
сохранить screen transforms для проверки фиксированных axes.

## energy_history.csv

```
time_s,kinetic_J_per_m,potential_J_per_m,capillary_J_per_m,wall_J_per_m,
physical_dissipation_integral_J_per_m,external_work_J_per_m,
numerical_energy_remainder_J_per_m,budget_defect_J_per_m,state_id
```

Если solver использует per-density units, преобразовать явно или назвать столбцы
соответствующе и записать units в manifest. Нельзя смешать разные соглашения.
При sigma=0 capillary term равен нулю по модели, не «измерен» как малая величина.
Если adsorption energy отсутствует как ведущепорядковое допущение, написать это;
не добавлять произвольное отрицательное тепло ради замыкания бюджета.

## runup_summary.json / refinement.json

Первый: есть ли новый левый максимум, его высота/время/неопределённость,
есть ли последующий более высокий, финальные H, какой record больше исходного P2.
Если события нет — null + reason, а не выдуманные нули.

Refinement: ID каждого run, единая physical config, изменённые mesh/dt/regularization,
общие времена, нормы различий surface/R/H, first peak metrics, cost, выводы.
Две точки уточнения подтверждают согласованность, но не доказывают асимптотический
порядок. Для claimed order нужны достаточные уровни и применимый режим.

## manifest.json

Минимально: schema_version, run_id, status, physical_model_label,
source_hash, execution_HEAD, resolved_case, artifacts[{path,sha256,role}],
native_data_paths, time_range, frame_map, contact_source_map, units, limitations.
Не включать собственный hash в циклическую самоссылку. Нельзя наследовать HEAD
исторического dataset как текущий execution HEAD.

## verification.json

Раздельно: contract, physics, wetting, numerical_refinement, rendering, artifact_integrity,
independent_review. У каждого status и evidence_paths. NOT_RUN не трактуется PASS.
Финальное COMPLETE разрешено только после полной обязательной проверки.

## Обязательные отрицательные тесты конечного verifier

1. H уменьшился на позднем кадре — отклонить.
2. H монотонен, но вырос без source reach — отклонить.
3. Старая P2 переместилась вверх/вниз — отклонить.
4. P3 исчезает после отступления волны — отклонить.
5. На W есть сухая дырка/разрыв слоя — отклонить.
6. Прямая в t=0 подменена кривой — отклонить.
7. Fixed endpoint + crest extrapolation помечены как actual contact — отклонить.
8. Смешаны run IDs/physical config между данными и refinement — отклонить.
9. Нет конечных данных t=5, но есть длинное циклическое MP4 — отклонить.
10. H пересчитан только по frame sampling с потерей внутреннего пика — отклонить.
11. Сохранённый dry/wet state расходится с нарисованной полосой — отклонить.
12. Синтетический fixture выдан за scientific output — отклонить.

Доступный в пакете tools/check_wetting_csv.py покрывает только часть пунктов
1–4 и геометрическую согласованность; остальные реализует агент в target verifier.
