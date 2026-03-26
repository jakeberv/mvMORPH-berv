# Specs

This folder holds durable project context for the local `mvMORPH` fork.

The goal is to make future threads fast to start and hard to derail:

- keep one short index of what exists
- keep one reusable template
- keep one detailed "fork delta" spec describing what changed relative to upstream
- keep one current-state spec describing what is validated, what is experimental, and what remains open

## Recommended Reading Order

For a new thread, read these first:

1. `specs/001-repo-context.md`
2. `specs/010-fork-delta.md`
3. `specs/020-validation-and-open-questions.md`
4. `specs/030-bmm-simulation-runbook.md` or `specs/040-oum-simulation-runbook.md`, depending on the workstream

## Conventions

- Prefer numbered files so reading order is obvious.
- Keep specs practical and current rather than exhaustive.
- When behavior changes, update the relevant spec in the same branch if possible.
- Treat specs as branch-local truth for the current fork, not as upstream package documentation.
- Prefer linking to concrete files and reports rather than pasting large logs.

## Current Spec Set

- `specs/001-repo-context.md`
  - repo layout, remotes, important directories, and current branch/worktree assumptions
- `specs/010-fork-delta.md`
  - detailed summary of all major changes built on top of the upstream fork baseline
- `specs/020-validation-and-open-questions.md`
  - what has been validated locally and remotely, plus the main unresolved questions
- `specs/030-bmm-simulation-runbook.md`
  - how to resume local or remote simulation work for the experimental corrpower BMM family
- `specs/040-oum-simulation-runbook.md`
  - how to resume local or remote simulation work for the OUM extensions and model-selection studies
- `specs/TEMPLATE.md`
  - starting point for future specs

## Scope

These specs are intentionally focused on the work added in this fork:

- `mvgls` OUM covariate support and associated validation
- experimental BMM correlation-aware work, culminating in `corrpower`
- simulation harnesses and summary reports added to support those changes

They are not a replacement for:

- the package `README`
- Rd help pages under `man/`
- simulation result reports under `reports/` and `tests/reports/`
