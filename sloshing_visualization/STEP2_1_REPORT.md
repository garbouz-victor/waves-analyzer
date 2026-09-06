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

Разделы с измерениями, исправленной геометрией и результатами проверок будут
добавлены после фактического postprocessing/rendering; здесь не заявлен их успех заранее.
