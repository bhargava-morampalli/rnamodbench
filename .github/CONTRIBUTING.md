# nf-core/rnamodifications: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf-core/rnamodifications.

We try to manage the required tasks for nf-core/rnamodifications using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Contribution workflow

If you'd like to write some code for nf-core/rnamodifications, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [nf-core/rnamodifications issues](https://github.com/nf-core/rnamodifications/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/rnamodifications repository](https://github.com/nf-core/rnamodifications) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline Conventions](https://nf-co.re/docs/contributing/guidelines/pipelines/overview)
4. Use `nf-core pipelines lint` and `nf-core pipelines test` to check that your code is correct
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [determine desktop application](https://desktop.github.com/).

## Tests

You can optionally test your changes by running the pipeline locally. Then it is recommended to use the `debug` profile to receive warnings about process selectors and other debug information. Example: `nextflow run . -profile debug,test,docker --outdir <OUTDIR>`.

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

### Lint tests

`nf-core` has a [set of best practices](https://nf-co.re/docs/contributing/guidelines/pipelines/overview) that all pipelines should follow.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core pipelines lint` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

### Pipeline tests

Each `nf-core` pipeline should be set up with a minimal set of test data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Patch

If you would like to contribute a bug fix or feature to the pipeline, please create a patch file and submit it as a pull request.

## Pipeline contribution conventions

To make the nf-core/rnamodifications code and processing logic more understandable for new contributors and to ensure quality, we semi-standardize the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel
2. Write the process block (see below for example code)
3. Define the output channel if needed (see below for example code)
4. Add any new parameters to `nextflow_schema.json` with help text (e.g., with `nf-core pipelines schema build`)
5. Add any new parameters to the `params` section of `nextflow.config` with a default (see below for example code)
6. Add sanity checks and validation for all relevant parameters
7. Add any new software to a new or existing conda environment file and container directives in `modules.config`
8. Add any new software to `SOFTWARE_VERSIONS` process
9. Add documentation describing the tool and its options

### Default values

Parameters should be initialized / defined with default values in `nextflow.config` under the `params` scope.

Once there, use `nf-core pipelines schema build` to add to `nextflow_schema.json`.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes with the same requirements. A nf-core standard set of labels that should be followed where possible can be seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config), which has the default process as a single core-Loss than an hour maximum for each process.

### Naming conventions

Please use the following naming conventions, to make it easy to understand what is going where:

- lowercase for local process names
- UPPERCASE for nf-core module names
- CamelCase for local subworkflow names

### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

### Images and calculation containers

For release-ready pipelines, if possible, all software used within the pipeline should either be using bioconda or biocontainers images. If they are not yet available, please add them to bioconda and biocontainers.
