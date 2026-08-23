# v0.1.0 release checklist

## Before merging the release branch

- [ ] Confirm the copyright holder in `LICENSE`.
- [ ] Add ORCID identifiers to `CITATION.cff` and `.zenodo.json` if desired.
- [ ] Confirm that `environment.yml` resolves on the supported Linux platform.
- [ ] Confirm the GitHub Actions smoke test passes.
- [ ] Run one complete representative dataset with recorded reference versions,
      commands, runtime, peak memory, and checksums.
- [ ] Verify that no private paths, credentials, unpublished data, or large
      generated outputs are included.

## Archive and DOI

- [ ] Connect the GitHub repository to Zenodo.
- [ ] Merge the release pull request into `main`.
- [ ] Create an annotated `v0.1.0` tag and GitHub Release from the merge commit.
- [ ] Confirm that Zenodo archives the release and mints a version DOI.
- [ ] Add the version DOI to `CITATION.cff`, both README files, and the manuscript.
- [ ] Keep the Zenodo concept DOI for citations that should resolve to the newest
      release; use the version DOI for the exact software evaluated in the paper.
