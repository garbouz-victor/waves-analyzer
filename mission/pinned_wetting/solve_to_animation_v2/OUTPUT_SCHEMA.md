# Минимальные экспортные интерфейсы внешних guards

Это адаптер к native, не новый обязательный внутренний storage format.
Все числа SI; `time_s` — физическое время. CSV comma, UTF-8 header, finite floats.
Дополнительные столбцы разрешены. Схемы рассчитаны на связную границу с одним
текущим контактом на каждой стенке; новая топология требует расширения **явно**,
не потери данных. Детектор обязан обрабатывать полный новый контакт с жидкой
областью, а не pixels. Область применения схемы не является доказательством,
что топология обязана оставаться такой 5 с.

## resolved_case.json

    {
      "physical_contract": <точная копия CONTRACT.json.physical_contract>,
      "numerics": {"method": "...", "source_hash": "...", "mesh": "..."},
      "provenance": {"actual_execution_HEAD": "..."}
    }

Guard проверяет declared fields. Фактические BC проверяет независимый native audit.

## contacts.csv — каждый accepted state, не render frames

    step,time_s,R_left_m,R_right_m,H_left_m,H_right_m,right_u_max_m_s,right_w_max_m_s,bottom_speed_max_m_s,left_normal_max_m_s

step 0,1,...; time строго возрастает от0 до5. R_right — правый **материальный**
контакт полной разрешённой границы, не уровень bulk на отступе от стены.
H=running max accepted endpoint samples; межшаговый экстремум поддерживается
необязательными `reach_left_m,reach_right_m` с отдельной native evidence.
Для t0 reach=R. В этом контракте справа R=H=z2. Нулевые BC поля — измеренные
trace max, не константы, написанные экспортёром по названию BC.

## retained_intervals.csv — обе стены в каждом состоянии

    step,side,z_bottom_m,z_top_m

side L/R; ровно один интервал [-d,H_i] каждой стороны каждого accepted step.

## interface_samples.csv — все 11 временных срезов, допустимы дополнительные

    time_s,vertex_id,x_m,z_m

Для каждого времени упорядоченные vertex_id=0..N-1; полная кривая от left контакта
к right P2. x может быть немонотонным: нельзя удалять почти вертикальный/обратный
участок. Это контролируемая геометрическая дискретизация native curve, не линейная
интерполяция картинки. Для криволинейного native сообщить tessellation error.

CSV snapshots должны совпадать с accepted временем или корректно привязанной
проверенной интерполяцией. Дешёвый checker пока требует совпадения со временем
contacts; для dense output экспортёр добавляет допустимый state с явным label
и не выдаёт его за accepted PDE step. Предпочтительно попадать dt в output times.

## marker_catalog.csv

    marker_id,side,x_m,z_m,created_at_s,source_step

P1/P2 обязательны, source_step=0,created_at_s=0. Все координаты fixed. Для новых
left рекордов source_step связывает marker height с реальным R того состояния;
created_at не раньше source time. Guard не требует P3/P4 и не доказывает, что
данный sample — локальный пик; это задача native detector review. Новых right
рекордов при фиксированной P2 в данной топологии нет.

## Проверка

    python mission/pinned_wetting/solve_to_animation_v2/tools/contract_guard.py OUTPUT/resolved_case.json
    python mission/pinned_wetting/solve_to_animation_v2/tools/check_exports.py OUTPUT

Exit0 означает **только** соответствие declarations/exports; CFD success по нему
не объявлять. Tests содержат synthetic fixtures, явно не physical data.
