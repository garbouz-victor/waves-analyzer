# Минимальная обвязка выполнения

## Использовать существующее

В текущем пакете PW1 уже есть CLI, lock, бюджет, checkpoints, renderer и review.
Новый scheduler/агентный framework не нужен. Добавить только новый model mode,
адаптер native/export для его геометрии и проверки. Каталог этой миссии хранит
постоянный контракт, текущие progress/решения живут в старом EXECPLAN/state с
явным новым scope или в одном связанном successor; никаких двух активных ledger.

Не изменять master, исторические native, STOP, symmetric output, старые contracts.
Рабочая ветка feature/pinned-wetting-history сейчас ac418484 (перепроверить перед
run); base commit — provenance, не обязательный текущий HEAD во время resume.
Новый contract приложением уточняет прежние инструкции только для нового mode.

## Ресурсы

Наследовать существующий общий предел 12 ч jobs/40 GiB и реальную уже списанную
сумму в state.json. Не начинать новые 12 ч заново. Перед запуском вычислить остаток.
Один heavy job, lock, PID/process-start-identity, heartbeat и checkpoint.
Платные облака, external API compute, auto-push, force/reset запрещены.
Сессия Codex, PDE process и сохранённый CFD checkpoint — разные уровни восстановления.
Не запускать второй solver, если процесс первого существует. Не обещать фоновой
работы без реально запущенного и зарегистрированного процесса.

Оценка времени рекомендательная, actual cap обязателен. Внутреннее время agent
reasoning и время solver jobs не смешивать; честно учитывать оба доступных счётчика.
Перед ресурсным STOP сохранить данные. Не ждать отсутствующий платный API.

## Проверяемые этапы

1. Контракт/среда/метод.
2. Новый solver, минимальные аналитические/геометрические тесты.
3. Pilot и preview с настоящей движущейся границей.
4. Target 5s, dt/h comparisons.
5. Renderer, native audit, release.

После PASS переход дальше. После solver failure — гипотеза/целевая правка,
не новый ручной запрос. До двух различающихся исправлений blocker разумно
выполнить в этой же миссии, если хватает ресурса; адаптивные шаги имеют свою
конечную retry policy. Не превращать лимит попыток в обязанность повторить
заведомо опасный job. Все отклонённые candidates сохраняются отдельно.

## Как не смешивать научную и административную identity

Численная identity: actual equations, BC per wall, физические параметры,
initial geometry, discretization/time controls, relevant source and library versions.
Execution lineage: actual HEAD, command, time, session, host/container.
Совместимое продолжение допускает другой HEAD при прежних numerical bytes.
Несовместимое численное изменение — новый run, не дописывание в старый native.

## Дешёвые guards пакета

`tools/contract_guard.py resolved_case.json`: сравнивает объявленные физические
поля с CONTRACT.json. Пример структуры создаётся через --example. Не читает solver.
`tools/check_exports.py OUTPUT`: проверяет CSV инварианты начальной прямой,
правой неподвижной P2, history H/W и каталога маркеров. Это НЕ PDE verifier.
Они stdlib-only, доступны и в чистом Python. Codex адаптирует exporter, не guard
под желаемый результат; schema в OUTPUT_SCHEMA.md.

## Git и установка

install.py внешнего пакета только копирует этот новый каталог, никогда не
переписывает отличающийся файл. --link-agents добавляет короткую ссылку;
изменение root AGENTS видно пользователю. Нет switch/reset/commit/push/solve.
После первого local commit повторить doctor/compatible restart regression.
GitHub CI утверждать только при реально наблюдаемом запуске; отсутствие push
значит, что локальные тесты не названы remote PASS.
