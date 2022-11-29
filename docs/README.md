# utilities docs

## Requirements

Install [`quarto`](https://quarto.org/docs/get-started/).

## Preview

Run `quarto preview docs` from the root of the repository to preview the website.

## Render

Run `quarto render docs` from the root of the repository to render the website.

## Module docs
The directory `modules/` was generated using:

```bash
bin/tools/native/quarto/generate_documentation_qmd/generate_documentation_qmd \
  --input "`pwd`/src" \
  --output "`pwd`/docs/modules" \
  --clean \
  --write_index \
  --git_repo "czbiohub/utilities" \
  --git_tag "main_build" \
  --git_browse_url "https://github.com/czbiohub/utilities/blob/main/"

bin/tools/docker/quarto/generate_documentation_qmd/generate_documentation_qmd \
  --input "`pwd`/module_openpipeline" \
  --output "`pwd`/docs/modules" \
  --git_repo "openpipelines-bio/openpipeline" \
  --git_tag "main_build" \
  --git_browse_url "https://github.com/openpipelines-bio/openpipeline/blob/main/"
```