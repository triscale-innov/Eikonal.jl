# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed - BREAKING CHANGES

- Removed `from_png` and convenience `FastSweeping` & `FastMarching`
  constructors from images. The same feature is now available via
  `Eikonal.img2array` (#22)


## [0.1.1] - 2023-08-21

### Changed

- Improve precompilation with PrecompileTools and a custom workload (#19)
- Improve type stability (#20)
- Update compat bounds (#21)


## [0.1.0] - 2023-08-17

First released version

### Added

- Solver based on the Fast Marching Method (FMM) for 2D Eikonal equations
- Solver based on the Fast Sweeping Method (FSM) for Eikonal equations in
  arbitrary dimension (#13)