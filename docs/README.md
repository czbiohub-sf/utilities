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
bin/tools/docker/quarto/generate_documentation_qmd/generate_documentation_qmd \
  --input "src" \
  --output "docs/modules" \
  --git_repo "czbiohub/utilities" \
  --git_tag "main_build" \
  --git_browse_url "https://github.com/czbiohub/utilities/blob/main/"
```