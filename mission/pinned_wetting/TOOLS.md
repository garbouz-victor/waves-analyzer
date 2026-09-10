# Что уже исполняется в этом пакете

## Проверка комплекта

```bash
python3 mission/pinned_wetting/tools/check_pack.py
```

Возвращает PACKAGE_OK для наличия/синтаксиса документов и скриптов. При наличии
PyYAML проверяет YAML и согласованность ключевых полей. Без PyYAML YAML-check
помечается SKIPPED. Ни одна ветвь не заявляет, что жидкость уже посчитана.

## Проверка wetting CSV

После того как агент создаст реальные файлы из OUTPUT_SCHEMA:

```bash
python3 mission/pinned_wetting/tools/check_wetting_csv.py \
  sloshing_visualization/output/pinned_wetting \
  --z1=-0.03492076949174773 --z2=0.03492076949174773 \
  --a=1 --t-end=5
```

Числа выше относятся **только к примеру** a=1, alpha=2°, z_mean=0.
В реальном run брать z1/z2 из resolved_case.json. Скрипт проверяет:
положительный horizon, initial state, running maximum с interval reach,
неизменные маркеры и их постоянное присутствие на кадрах, связи state IDs.

Он не проверяет PDE, реальную связь R с полями, вычисление плёнки, convergence
или соответствие MP4 tracks. Вывод явно содержит NOT_VALIDATED_BY_THIS_CHECK.
Агент должен реализовать полный scientific/visual verifier поверх этого subset.

## Unit tests

```bash
python3 -m unittest discover -s mission/pinned_wetting/tests -v
```

Все тестовые ряды SYNTHETIC и создаются только во временной папке. Их нельзя
использовать для final_animation. В пакете нет example scientific results.

## Изменение checker

Если из-за обоснованного метода требуется расширение CSV, сохранить старые поля
или сделать явную новую версию схемы. Не удалять отрицательные тесты, чтобы получить
PASS. Source contact heights всё равно независимо сверяются по native solver data.
