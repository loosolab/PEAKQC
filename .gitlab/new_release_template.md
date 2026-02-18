# Checklist for new version release

## Before merge:
- [ ] Update date in Changes.md
- [ ] Update peakqc/_version.py

## After merge:
### Gitlab (pipeline takes ~1.5 hrs)
- [ ] Create new release (automated)
- [ ] Add new version tag (automated)
- [ ] Add milestone to release and close it (if one is available)
### Github (https://github.com/loosolab/PEAKQC)
The main and dev branch are automatically mirrored.
- [ ] check if the repository updated (may take a few minutes)
- [ ] Create a new release (copy from Gitlab)
### PyPI
The final step in the CI/CD pipeline creates a new PyPI release (after ~1.5 hrs).
- [ ] check if the release was sucessful
