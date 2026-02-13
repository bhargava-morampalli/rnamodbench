# bhargava-morampalli/rnamodbench: Contributing Guidelines

Thanks for contributing to `rnamodbench`.

## Workflow

1. Check existing issues before starting work.
2. Fork the repository and create a branch from `dev`.
3. Implement your changes with tests when relevant.
4. Run local checks before opening a PR:
   - `nextflow config`
   - `nf-test test` (or targeted nf-tests)
   - `pytest` for Python components when applicable
5. Open a pull request against `dev` with a clear summary.

## Expectations

- Keep changes scoped and documented.
- Update `docs/usage.md` and `docs/output.md` if behavior changes.
- Update `CHANGELOG.md` for user-facing changes.
- Avoid introducing hard-coded machine-specific paths.

## Pipeline conventions

- Use lowercase for local process names.
- Use uppercase for module process identifiers.
- Keep resource labels aligned with `conf/base.config`.
- Add new params to `nextflow_schema.json` and `nextflow.config`.

## Need help?

Open a discussion or issue in this repository:
`https://github.com/bhargava-morampalli/rnamodbench`.
