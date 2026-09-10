# Происхождение требований и источники

Подготовлено 2026-09-10. Репозиторий прочитан через GitHub connection, а не по
предположению о содержимом файлов. Ни один репозиторный файл этим пакетом не изменён.

## Пользовательские требования

Основной источник — последние сообщения текущего разговора, процитированные
в PHYSICS.md: необратимое присутствие жидкости и непрерывный остаточный слой
ниже максимума; точная толщина не нужна; начальная прямая; фиксированные координаты
точек; возможность повысить точность расчётом. a=1, d=10, срезы 0–5 с через 0.5 с
наследуются из более ранней исходной задачи.

Nu/alpha не были однозначно выбраны пользователем для новой миссии. Значения
0.01 м²/с и 2° — авторский стартовый пример из defaults проекта. Бюджет 12 часов
и 40 GiB — организационное предложение пакета. Нет гарантии завершения за этот срок.

Раньше ассистент предлагал экстраполяцию near-wall fixed-endpoint поля как run-up.
В этом пакете она **не допускается как физически квалифицированный конечный путь**.

## Проверенные файлы проекта

Repository: https://github.com/garbouz-victor/waves-analyzer
Observed master: 7f0c45f192520b218838d535b0182facf43133d3.

1. Linear config, включая default parameters:
https://github.com/garbouz-victor/waves-analyzer/blob/7f0c45f192520b218838d535b0182facf43133d3/sloshing_visualization/src/sloshing/config.py

2. FEM trace/no-slip endpoint structure:
https://github.com/garbouz-victor/waves-analyzer/blob/7f0c45f192520b218838d535b0182facf43133d3/sloshing_visualization/src/sloshing/fem_spaces.py

3. Surface time stepping:
https://github.com/garbouz-victor/waves-analyzer/blob/7f0c45f192520b218838d535b0182facf43133d3/sloshing_visualization/src/sloshing/time_integrator.py

4. Full rate restriction: matched density, zero forcing:
https://github.com/garbouz-victor/waves-analyzer/blob/7f0c45f192520b218838d535b0182facf43133d3/sloshing_visualization/src/sloshing/multiphase/full_phase_rate.py

5. Rendering and declared limitations:
https://github.com/garbouz-victor/waves-analyzer/blob/7f0c45f192520b218838d535b0182facf43133d3/sloshing_visualization/README.md

## Внешний технический контекст — не дополнительные пользовательские требования

Официальные материалы OpenAI о постоянных instructions, execution plans и
возобновлении CLI. Использованы для организации задания, не как доказательство CFD:

https://developers.openai.com/codex/guides/agents-md/
https://developers.openai.com/cookbook/articles/codex_exec_plans/
https://developers.openai.com/codex/noninteractive/

Установленный codex --help является источником истины для доступных флагов.
Пакет не требует особой версии CLI, платного API key или стороннего orchestration SDK.
Не выдавать capability/session limits за обещание бессрочной автономной работы.

Физический контекст: D. N. Sibley, A. Nold, N. Savva, S. Kalliadasis,
“The contact line behaviour of solid-liquid-gas diffuse-interface models”:
https://arxiv.org/abs/1310.1255

Работа иллюстрирует различие sharp/diffuse descriptions с no-slip. Она не доказывает
нашу модель «жидкость приклеивается навсегда» и не даёт недостающие параметры стены.
Эта необратимость — конкретное допущение пользователя, нуждающееся в собственной
численной реализации и проверке.
