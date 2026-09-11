# Источники и пределы их применимости

Проверено 2026-09-10. Это references для работы агента, не утверждение, что точная
задача уже решена в какой-либо из статей. В следующем исследовании использовать
полный текст релевантных первоисточников; из одного abstract не выводить методы.

## Код пользователя: именно источник текущих параметров/ограничений

- https://github.com/garbouz-victor/waves-analyzer/blob/ac418484d48dff5f8104d2177b5c2b1f703a232a/mission/pinned_wetting/left_navier_right_noslip/model_decision.md
- https://github.com/garbouz-victor/waves-analyzer/blob/ac418484d48dff5f8104d2177b5c2b1f703a232a/sloshing_visualization/src/sloshing/pinned_wetting/asymmetric.py
- https://github.com/garbouz-victor/waves-analyzer/blob/ac418484d48dff5f8104d2177b5c2b1f703a232a/mission/pinned_wetting/left_navier_right_noslip/review.md

Здесь fixed-domain linear physics, неподвижная P2 и провал all-edge slope,
не nonlinear solver. Новое разрешение full NS явно в MASTER_PROMPT.

## Кинематика и численные методы

- Fricke, Köhne, Bothe, A Kinematic Evolution Equation for the Dynamic Contact
  Angle and some Consequences, Physica D 394 (2019),26–43.
  https://arxiv.org/abs/1810.00830 ; DOI 10.1016/j.physd.2019.01.008.
  Связь contact-angle dynamics с transporting velocity/regularity. Не является
  готовым решением сосуда пользователя; Navier не гарантирует smooth corner.
- Fricke, Marić, Bothe, Contact line advection using the geometrical VOF method,
  JCP407 (2020),109221. https://arxiv.org/abs/1907.01785 .
  Проверка boundary interface transport важна сама по себе. Prescribed velocity
  advection benchmark не заменяет coupled momentum/gravity solve.
- Guo,Tice, Decay of viscous surface waves without surface tension.
  https://arxiv.org/abs/1011.5179 . Free moving surface с sigma=0 — содержательная
  модель; геометрии статьи periodic/infinite, не доказательство нашего бокового
  pinned-contact случая и не гарантия его глобальной регулярности.
- Hu,Lei,Li,Tang, Energy Dissipating ALE-MDR Method for Navier–Stokes Free Boundary
  Problems with Moving Contact Line (2026), DOI10.1137/25M1784958.
  https://epubs.siam.org/doi/10.1137/25M1784958 . Идеи ALE/mesh quality; статья
  содержит capillary/prescribed-angle physics: не переносить их в sigma=0 target.

## Организация автономного Codex

- OpenAI, Using PLANS.md for multi-hour problem solving:
  https://developers.openai.com/cookbook/articles/codex_exec_plans
  Living execution plan и переход через milestones без повторных запросов.
- OpenAI, Harness engineering:
  https://openai.com/index/harness-engineering/
  Короткий AGENTS как карта; результат проверяется средой, не длинной декларацией.
- Official Codex instructions:
  https://developers.openai.com/codex/guides/agents-md/
  https://developers.openai.com/codex/noninteractive/
  Конкретный запуск/флаги сверять с installed `codex --help`; пакет не устанавливает
  Codex и не запускает бесконечный loop. Активный user prompt явно называет новую
  миссию, старые линейные instructions остаются historical.
